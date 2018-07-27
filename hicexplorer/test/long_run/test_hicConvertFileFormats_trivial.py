from hicexplorer import hicConvertFileFormats

def test_help():
    args = "--help"
    hicConvertFileFormats.main(args)

def test_version():
    args = "--version"
    hicConvertFileFormats.main(args)