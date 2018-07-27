from hicexplorer import hicSumMatrices

def test_help():
    args = "--help"
    hicSumMatrices.main(args)

def test_version():
    args = "--version"
    hicSumMatrices.main(args)