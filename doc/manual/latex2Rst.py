import re
import sys
import os

try:
    filename = sys.argv[1]
    f = open(filename, 'r')
    textIn = f.read()
    f.close()
    filename = filename.split('.')[0]+".rst"
except:
    filename = 'output.rst'
    print ">>> The script will use default test: search and replace the inline math mode for the following text"

    # ...
    textIn = 'in this section, $P=(x,y, z)$ will denote the coordinate of the control point'
    textIn += ' '
    textIn += 'the solutions of this equation are given by \cite{Fermat1680} blablabla'
    textIn += ' '
    textIn += 'to see the effect of moving Ctrl Pts see Figure \\ref{fig:move_ctrl_pts}'
    # ...

def inlineEnvMath2(s):
    """
    replace the inline latex math mode $ ... $ by the Sphinx one :math:`...`
    """
    rule = '\$\$(.*?)\$\$'
    leftLim = ':math:`'
    rightLim = '`'

    p = re.compile(rule)
    #p = re.compile(rule, re.IGNORECASE)

    # ...
    list_data = p.split(s) # split the whole text with respect to the rule
    list_exp  = p.findall(s) # get all expressions to replace

    list_text = ""
    for data in list_data:
        new_data = data
        if data in list_exp:
            new_data = leftLim + data + rightLim
        list_text += new_data

    return list_text

def inlineEnvMath(s):
    """
    replace the inline latex math mode $ ... $ by the Sphinx one :math:`...`
    """
    rule = '\$(.*?)\$'
    leftLim = ':math:`'
    rightLim = '`'

    p = re.compile(rule)
    #p = re.compile(rule, re.IGNORECASE)

    # ...
    list_data = p.split(s) # split the whole text with respect to the rule
    list_exp  = p.findall(s) # get all expressions to replace

    list_text = ""
    for data in list_data:
        new_data = data
        if data in list_exp:
            new_data = leftLim + data + rightLim
        list_text += new_data

    return list_text

#def equationEnvMath(s):
#    """
#    TODO : Not yet validated
#    replace the inline latex math mode \\begin{equation} ... \end{equation} by the Sphinx one :math:`...`
#    """
#    TAG = 'equation'
#    beginLim = '\\begin{' + TAG + '}'
#    endLim   = '\end{'    + TAG + '}'
#    rule = '\\' + beginLim + '(.*?)' + '\\' + endLim
#    leftLim = '<<'
#    rightLim = '>>'
#
#    p = re.compile(rule)
#    #p = re.compile(rule, re.IGNORECASE)
#
#    # ...
#    list_data = p.split(s) # split the whole text with respect to the rule
#    list_exp  = p.findall(s) # get all expressions to replace
#
#    list_text = ""
#    for data in list_data:
#        new_data = data
#        if data in list_exp:
#            new_data = leftLim + data + rightLim
#        list_text += new_data
#
#    return list_text

def inlineEnvRef(s):
    """
    replace the inline latex ref mode \ref{...} by the Sphinx one :ref:`...`
    """
    rule = '\\\ref{(.*?)}'
    leftLim = ':ref:`'
    rightLim = '`'

    p = re.compile(rule)
    #p = re.compile(rule, re.IGNORECASE)

    # ...
    list_data = p.split(s) # split the whole text with respect to the rule
    list_exp  = p.findall(s) # get all expressions to replace

    list_text = ""
    for data in list_data:
        new_data = data
        if data in list_exp:
            new_data = leftLim + data + rightLim
        list_text += new_data

    return list_text

def inlineEnvCite(s):
    """
    replace the inline latex cite mode \cite{...} by the Sphinx one :cite:`...`
    """
    rule = '\\\cite{(.*?)}'
    leftLim = ':cite:`'
    rightLim = '`'

    p = re.compile(rule)
    #p = re.compile(rule, re.IGNORECASE)

    # ...
    list_data = p.split(s) # split the whole text with respect to the rule
    list_exp  = p.findall(s) # get all expressions to replace

    list_text = ""
    for data in list_data:
        new_data = data
        if data in list_exp:
            new_data = leftLim + data + rightLim
        list_text += new_data

    return list_text

def inlineEnvCitep(s):
    """
    replace the inline latex cite mode \citep{...} by the Sphinx one :cite:`...`
    """
    rule = '\\\citep{(.*?)}'
    leftLim = ':cite:`'
    rightLim = '`'

    p = re.compile(rule)
    #p = re.compile(rule, re.IGNORECASE)

    # ...
    list_data = p.split(s) # split the whole text with respect to the rule
    list_exp  = p.findall(s) # get all expressions to replace

    list_text = ""
    for data in list_data:
        new_data = data
        if data in list_exp:
            new_data = leftLim + data + rightLim
        list_text += new_data

    return list_text

def inlineEnvTextit(s):
    """
    replace the inline latex cite mode \textit{...} by the Sphinx one
    """
    rule = '\\\textit{(.*?)}'
    leftLim  = '*'
    rightLim = '*'

    p = re.compile(rule)
    #p = re.compile(rule, re.IGNORECASE)

    # ...
    list_data = p.split(s) # split the whole text with respect to the rule
    list_exp  = p.findall(s) # get all expressions to replace

    list_text = ""
    for data in list_data:
        new_data = data
        if data in list_exp:
            new_data = leftLim + data + rightLim
        list_text += new_data

    return list_text


def inlineEnvTextbf(s):
    """
    replace the inline latex cite mode \textbf{...} by the Sphinx one
    """
    rule = '\\\textbf{(.*?)}'
    leftLim  = '**'
    rightLim = '**'

    p = re.compile(rule)
    #p = re.compile(rule, re.IGNORECASE)

    # ...
    list_data = p.split(s) # split the whole text with respect to the rule
    list_exp  = p.findall(s) # get all expressions to replace

    list_text = ""
    for data in list_data:
        new_data = data
        if data in list_exp:
            new_data = leftLim + data + rightLim
        list_text += new_data

    return list_text

def inlineEnvTexttt(s):
    """
    replace the inline latex cite mode \texttt{...} by the Sphinx one :cite:`...`
    """
    rule = '\\\texttt{(.*?)}'
    leftLim  = '``'
    rightLim = '``'

    p = re.compile(rule)
    #p = re.compile(rule, re.IGNORECASE)

    # ...
    list_data = p.split(s) # split the whole text with respect to the rule
    list_exp  = p.findall(s) # get all expressions to replace

    list_text = ""
    for data in list_data:
        new_data = data
        if data in list_exp:
            new_data = leftLim + data + rightLim
        list_text += new_data

    return list_text

def inlineEnvFigRef(s):
    """
    replace the inline latex cite mode \ref{fig:...} by the Sphinx one :num:`Fig. #fig-...`
    """
    rule = '\\\\ref{fig:(.*?)}'
    leftLim = ':num:`Fig. #fig-'
    rightLim = '`'

    p = re.compile(rule)
    #p = re.compile(rule, re.IGNORECASE)

    # ...
    list_data = p.split(s) # split the whole text with respect to the rule
    list_exp  = p.findall(s) # get all expressions to replace

    list_text = ""
    for data in list_data:
        new_data = data
        if data in list_exp:
            new_data = leftLim + data + rightLim
        list_text += new_data

    return list_text

def inlineEnvEqRef(s):
    """
    replace the inline latex cite mode \eqref{eq:...} by the Sphinx one :eq:`eq-...`
    """
    rule = '\\\eqref{eq:(.*?)}'
    leftLim = ':eq:`eq-'
    rightLim = '`'

    p = re.compile(rule)
    #p = re.compile(rule, re.IGNORECASE)

    # ...
    list_data = p.split(s) # split the whole text with respect to the rule
    list_exp  = p.findall(s) # get all expressions to replace

    list_text = ""
    for data in list_data:
        new_data = data
        if data in list_exp:
            new_data = leftLim + data + rightLim
        list_text += new_data

    return list_text

# ... validated
textIn = inlineEnvMath2(textIn)
textIn = inlineEnvMath(textIn)
textIn = inlineEnvRef(textIn)
textIn = inlineEnvCite(textIn)
textIn = inlineEnvCitep(textIn)
textIn = inlineEnvTexttt(textIn)
textIn = inlineEnvFigRef(textIn)
textIn = inlineEnvEqRef(textIn)
# ...

# ... not yet validated
#textOut = equationEnvMath(textIn)
# ...

f = open(filename, 'w')
f.write(textIn)
f.close()
# ...




