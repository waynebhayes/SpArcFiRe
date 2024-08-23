def vprint(verbosity, *args, **kwargs):
    if verbosity:
        print(*args, **kwargs)