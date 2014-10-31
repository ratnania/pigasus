# -*- coding: UTF-8 -*-
#! /usr/bin/python

__author__="ARA"
__all__ = ['pigasusObject']

from . import common_obj as com

class pigasusObject(object):
    def __init__(self, *args, **kwargs):
       self.com = com.common_obj()

       self.infoData = {} # dictionary for info data. used for printing

    def writeInfoDataXML(self, doc, rootElt):

        # Create the main <card> element
        CURRENT_TAG = self.__class__.__name__
        maincard = doc.createElement(CURRENT_TAG)
        rootElt.appendChild(maincard)

        for d in self.infoData:
            TAG = d ; txt = self.infoData[d]
            curElt = doc.createElement(TAG)
            curText = doc.createTextNode(txt)
            curElt.appendChild(curText)
            maincard.appendChild(curElt)

    def __str__(self):
        line = ""
        for d in self.infoData:
            line += str(d + " : " + str(self.infoData[d]))
            line += "\n"
        return line

    def __eq__(self, other):
        return self.id == other.id

    def __ne__(self, other):
        return self.id != other.id




