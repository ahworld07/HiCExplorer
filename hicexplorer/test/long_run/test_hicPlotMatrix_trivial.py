from hicexplorer import hicPlotMatrix

def test_help():
    args = "--help"
    hicPlotMatrix.main(args)

def test_version():
    args = "--version"
    hicPlotMatrix.main(args)