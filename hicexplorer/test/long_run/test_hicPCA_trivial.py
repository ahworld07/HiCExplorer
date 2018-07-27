from hicexplorer import hicPCA

def test_help():
    args = "--help"
    hicPCA.main(args)

def test_version():
    args = "--version"
    hicPCA.main(args)