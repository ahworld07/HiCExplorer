from hicexplorer import hicCorrectMatrix

def test_help():
    args = "--help"
    hicCorrectMatrix.main(args)

def test_version():
    args = "--version"
    hicCorrectMatrix.main(args)