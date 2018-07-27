from hicexplorer import hicFindEnrichedContacts


def test_help():
    args = "--help"
    hicFindEnrichedContacts.main(args)

def test_version():
    args = "--version"
    hicFindEnrichedContacts.main(args)