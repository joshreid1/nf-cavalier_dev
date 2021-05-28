
def read_tsv(Path path, List<String> names) {
    path.toFile().readLines().with { lines ->
        assert names == lines[0].split('\t')
        lines.drop(1).collect { [names, it.split('\t')].transpose()
            .collectEntries { k, v -> [(k): v] }
        }
    }
}