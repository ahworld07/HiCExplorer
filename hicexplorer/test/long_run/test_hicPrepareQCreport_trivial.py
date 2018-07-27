from hicexplorer import hicPrepareQCreport

def test_help():
    args = "--help"
    hicPrepareQCreport.main(args)

def test_version():
    args = "--version"
    hicPrepareQCreport.main(args)