from hicexplorer import hicPlotTADs

def test_help():
    args = "--help"
    hicPlotTADs.main(args)

def test_version():
    args = "--version"
    hicPlotTADs.main(args)