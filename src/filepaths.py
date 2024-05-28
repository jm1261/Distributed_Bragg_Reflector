import os


def check_dir_exists(dir_path : str) -> None:
    """
    Function Details
    ================
    Check to see if a directory path exists, and if not create one.

    Parameters
    ----------
    dir_path: string
        Path to target directory.

    Returns
    -------
    None.

    See Also
    --------
    os path isdir
    os mkdir

    Notes
    -----
    None.

    Example
    -------
    None.

    ----------------------------------------------------------------------------
    Update History
    ==============

    25/05/2024
    ----------
    Documentation updated.

    """
    if os.path.isdir(dir_path) is False:
        os.mkdir(dir_path)