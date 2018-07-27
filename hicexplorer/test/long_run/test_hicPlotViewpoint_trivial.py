from hicexplorer import hicPlotViewpoint

def test_help():
    args = "--help"
    hicPlotViewpoint.main(args)

def test_version():
    args = "--version"
    hicPlotViewpoint.main(args)