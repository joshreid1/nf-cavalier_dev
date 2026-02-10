
process PDF_SPLIT {
    label 'C2M2T2'
    label 'qpdf'
    tag "$fam"

    input:
    tuple val(fam), path(pdf)

    output:
    path("genes/*.pdf")

    script:
    """
    mkdir -p genes

    qpdf --json --json-key=outlines $pdf \\
      | awk '
        /"title":/ {
          ++p
          line=\$0

          if (line ~ /:[[:space:]]*[^[:space:]]+[[:space:]]*-[[:space:]]*chr/) {
            sub(/.*:[[:space:]]*/, "", line)
            sub(/[[:space:]]*-[[:space:]]*chr.*/, "", line)
            gene=toupper(line)

            if (!(gene in seen)) { seen[gene]=1; order[++n]=gene }

            # store pages in a 2D-ish array keyed by gene + index
            k = ++cnt[gene]
            page[gene, k] = p
          }
        }

        END {
          for (i=1; i<=n; i++) {
            gene = order[i]
            s = ""
            for (j=1; j<=cnt[gene]; j++) {
              if (j>1) s = s ","
              s = s page[gene, j]
            }
            print "qpdf --empty --pages \\"$pdf\\" " s " -- \\"" "genes/${fam}." gene ".pdf\\""
          }
        }
      ' \\
      | xargs -P ${task.cpus} -I{} bash -lc '{}'
    """
}