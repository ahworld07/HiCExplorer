from hicexplorer import hicInfo

def test_help():
    args = "--help"
    hicInfo.main(args)

def test_version():
    args = "--version"
    hicInfo.main(args)