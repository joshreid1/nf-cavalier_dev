
process pull_latest_list {
    label 'C1M1T1'
    label 'cavalier'
//    container = null
//    module = 'R/3.6.1'
    publishDir "output/pull_latest_list", mode: 'copy'
    tag { "$id:v$ver" }

    input:
        tuple val(id), val(ver)

    output:
        tuple val(id), path(output)

    script:
        output = id.replaceAll(':', '_') + '_v' + ver + '.tsv'
        cmd = "cavalier::get_web_list('$id', version='$ver', save='$output', secure=FALSE)"
        """
        R --slave --vanilla -e "$cmd"
        """
}

