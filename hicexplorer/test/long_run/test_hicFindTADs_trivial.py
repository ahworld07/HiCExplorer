from hicexplorer import hicFindTADs


def test_help():
    args = "--help"
    hicFindTADs.main(args)

def test_version():
    args = "--version"
    hicFindTADs.main(args)