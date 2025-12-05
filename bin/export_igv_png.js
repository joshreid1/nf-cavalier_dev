#!/usr/bin/env node

const fs = require("fs");
const path = require("path");
const puppeteer = require("puppeteer");

async function sleep(ms) {
    return new Promise(r => setTimeout(r, ms));
}

// ---------------------------------------------------------------------------
// NEW: robust wait-for-download helper
// ---------------------------------------------------------------------------
async function waitForNewPNG(dir, timeout = 15000) {
    const start = Date.now();
    const initial = new Set(
        fs.readdirSync(dir).filter(f => f.endsWith(".png"))
    );

    while (Date.now() - start < timeout) {
        const current = fs.readdirSync(dir).filter(f => f.endsWith(".png"));
        for (const f of current) {
            if (!initial.has(f)) {
                return path.join(dir, f);
            }
        }
        await sleep(100);
    }

    return null;  // timed out
}
// ---------------------------------------------------------------------------

(async () => {

    const outDir = ".";
    const input = process.argv[2];
    const prefix = process.argv[3] || "";

    const viewportWidth  = process.argv[4] ? parseInt(process.argv[4], 10) : 800;
    const viewportHeight = process.argv[5] ? parseInt(process.argv[5], 10) : 800;
    const scaleFactor    = process.argv[6] ? parseFloat(process.argv[6]) : 2;

    if (!input) {
        console.error("Usage: export_igv_png.js <report.html> [prefix] [width] [height] [scale]");
        process.exit(1);
    }

    fs.mkdirSync(outDir, { recursive: true });

    const browser = await puppeteer.launch({
        headless: "new",
        executablePath: process.env.CHROMIUM_PATH || undefined,
        defaultViewport: {
            width: viewportWidth,
            height: viewportHeight,
            deviceScaleFactor: scaleFactor
        },
        args: [
            "--no-sandbox",
            "--disable-setuid-sandbox",
            "--allow-file-access-from-files"
        ]
    });

    const page = await browser.newPage();

    // Enable downloads
    const client = await page.target().createCDPSession();
    await client.send("Page.setDownloadBehavior", {
        behavior: "allow",
        downloadPath: path.resolve(outDir)
    });

    console.log("Loading page:", input);
    await page.goto(`file:${path.resolve(input)}`, {
        waitUntil: "networkidle0"
    });

    // ---- WAIT FOR VARIANT TABLE ----
    console.log("Waiting for variant table…");
    await page.waitForFunction(() => {
        const table = document.querySelector("#variant_table");
        return table && table.querySelectorAll("tr").length > 1;
    }, { timeout: 60000 });

    console.log("Variant table detected.");

    // Extract variants list inside page
    const variants = await page.evaluate(() => {
        const rows = [...document.querySelectorAll("#variant_table tr")];
        return rows
            .map((row, idx) => {
                const tds = row.querySelectorAll("td");
                if (tds.length < 4) return null;
                return {
                    rowIndex: idx,
                    id: tds[3].innerText.trim()
                };
            })
            .filter(v => v && v.id !== "");
    });

    console.log(`Found ${variants.length} variants.`);

    // -----------------------------------------------------------------------
    //   FUNCTION: Trigger Save → Save as PNG (patched)
    // -----------------------------------------------------------------------
    async function savePNG(variantID) {

        await page.evaluate(() => {
            const root = document.querySelector("#igvDiv").shadowRoot;
            const btn = root.querySelector('.igv-navbar-icon-button[title="Save Image"]');
            ["pointerdown", "mousedown", "pointerup", "mouseup", "click"].forEach(ev =>
                btn.dispatchEvent(new Event(ev, { bubbles: true, composed: true }))
            );
        });

        await page.waitForFunction(() => {
            const root = document.querySelector("#igvDiv").shadowRoot;
            return !!root.querySelector(".igv-ui-dropdown");
        });

        await page.evaluate(() => {
            const root = document.querySelector("#igvDiv").shadowRoot;
            const pngItem = [...root.querySelectorAll(".igv-ui-dropdown div div")]
                .find(el => el.textContent.trim() === "Save as PNG");

            ["pointerdown", "mousedown", "pointerup", "mouseup", "click"].forEach(ev =>
                pngItem.dispatchEvent(new Event(ev, { bubbles: true, composed: true }))
            );
        });

        // -------------------------------------------------------------------
        // NEW: Wait until actual PNG appears (instead of fixed sleep)
        // -------------------------------------------------------------------
        const downloaded = await waitForNewPNG(outDir, 15000);
        if (!downloaded) {
            throw new Error(`PNG did not appear in ${outDir} after saving ${variantID}`);
        }

        const finalName = path.join(outDir, `${prefix}.${variantID}.png`);
        fs.renameSync(downloaded, finalName);
    }

    // -----------------------------------------------------------------------
    //   FUNCTION: click cell to update IGV locus
    // -----------------------------------------------------------------------
    async function clickVariantRow(idx) {
        await page.evaluate((rowIndex) => {
            const row = document.querySelectorAll("#variant_table tr")[rowIndex];
            if (!row) return;
            const cell = row.querySelectorAll("td")[3];
            if (!cell) return;

            ["pointerdown", "mousedown", "mouseup", "pointerup", "click"].forEach(ev =>
                cell.dispatchEvent(new Event(ev, { bubbles: true, composed: true }))
            );
        }, idx);
    }

    // -----------------------------------------------------------------------
    //   MAIN LOOP
    // -----------------------------------------------------------------------
    for (const v of variants) {
        console.log(`→ Selecting variant: ${v.id}`);

        await clickVariantRow(v.rowIndex);
        await sleep(100);   // give IGV time to render new region

        console.log("   Saving PNG…");
        await savePNG(v.id);

        console.log(`   ✔ Saved ${prefix}.${v.id}.png`);
    }

    console.log("All variants processed.");

    await browser.close();
})();
