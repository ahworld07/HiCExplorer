from hicexplorer import hicCorrelate

def test_help():
    args = "--help"
    hicCorrelate.main(args)

def test_version():
    args = "--version"
    hicCorrelate.main(args)