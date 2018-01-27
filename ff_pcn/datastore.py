#!/usr/bin/env python

"""
Module abstracting a datastore to hold information persitently.
"""


import os


class DataStore(object):

    def __init__(self):
        self.datastore = dict()

    def add(self, section, key, value):
        """
        Add key with value to section.
        Writes datatore file afterwards.
        """
        key = str(key)
        if not section in self.datastore:
            self.datastore[section] = dict()
        self.datastore[section][key] = value

    def get(self, section, key, default=None):
        """
        Exception save getting of key in section.
        """
        key = str(key)
        section = self.datastore.get(section, None)
        if section:
            return section.get(key, None)


datastore = DataStore()
"""Global datastore."""


def store(section):
    """
    Decorator to save returned values to datastore.
    """
    def store_decorator(func):
        def wrapped_func(*k):
            val = datastore.get(section, k)
            if val is not None:
                return val
            val = func(*k)
            datastore.add(section, k, val)
            return val
        return wrapped_func
    return store_decorator
