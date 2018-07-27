from hicexplorer import hicTransform

def test_help():
    args = "--help"
    hicTransform.main(args)

def test_version():
    args = "--version"
    hicTransform.main(args)