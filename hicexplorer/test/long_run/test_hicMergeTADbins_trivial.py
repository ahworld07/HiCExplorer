from hicexplorer import hicMergeTADbins

def test_help():
    args = "--help"
    hicMergeTADbins.main(args)

def test_version():
    args = "--version"
    hicMergeTADbins.main(args)