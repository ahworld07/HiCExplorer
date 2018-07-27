from hicexplorer import hicMergeMatrixBins


def test_help():
    args = "--help"
    hicMergeMatrixBins.main(args)

def test_version():
    args = "--version"
    hicMergeMatrixBins.main(args)