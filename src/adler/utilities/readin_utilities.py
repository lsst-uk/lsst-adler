import logging

logger = logging.getLogger(__name__)


def read_in_SSObjectID_file(readin_path):
    with open(readin_path) as f:
        ssoid_list = f.read().splitlines()

    # bit of validation here: we expect these to cast nicely to ints
    for ssoid in ssoid_list:
        try:
            int(ssoid)
        except ValueError:
            logger.error(
                "ValueError: One or more of the SSObjectIDs in the supplied list does not seem to be valid."
            )
            raise ValueError("One or more of the SSObjectIDs in the supplied list does not seem to be valid.")

    return ssoid_list
