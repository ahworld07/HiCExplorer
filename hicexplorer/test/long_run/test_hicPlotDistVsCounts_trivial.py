from hicexplorer import hicPlotDistVsCounts

def test_help():
    args = "--help"
    hicPlotDistVsCounts.main(args)

def test_version():
    args = "--version"
    hicPlotDistVsCounts.main(args)