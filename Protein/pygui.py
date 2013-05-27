from easygui import *

# currently can do one shape, one simulation at a time. will change later.
# add physical parameter info for each shape
# account for dummy values
# update protein.pl
# standardize file naming

shape = ['randst', 'sp', 'e', 'b', 'p', 'c', 'TIE_FIGHTER', 'triangle']
choice = choicebox(msg='Chosen shape:', title=' ', choices=shape)

if choice == 'randst':
    multenterbox(msg='Fill in values for the fields.', title=' ', fields=(), values=())
    msg = "Enter parameters:"
    title = "Simulation info."
    fieldNames = ["Length of cylinder:", "Endcap radius:"]
    fieldValues = []  # we start with blanks for the values
    fieldValues = multenterbox(msg,title, fieldNames)
    print "Reply was:", fieldValues
if choice == 'sp':
    multenterbox(msg='Fill in values for the fields.', title=' ', fields=(), values=())
    msg = "Enter your personal information"
    title = "Credit Card Application"
    fieldNames = ["Name","Street Address","City","State","ZipCode"]
    fieldValues = []  # we start with blanks for the values
    fieldValues = multenterbox(msg,title, fieldNames)
    print "Reply was:", fieldValues
if choice == 'TIE_fighter':
    multenterbox(msg='Fill in values for the fields.', title=' ', fields=(), values=())
    msg = "Enter your personal information"
    title = "Credit Card Application"
    fieldNames = ["Name","Street Address","City","State","ZipCode"]
    fieldValues = []  # we start with blanks for the values
    fieldValues = multenterbox(msg,title, fieldNames)
    print "Reply was:", fieldValues
if choice == 'triangle':
    multenterbox(msg='Fill in values for the fields.', title=' ', fields=(), values=())
    msg = "Enter your personal information"
    title = "Credit Card Application"
    fieldNames = ["Name","Street Address","City","State","ZipCode"]
    fieldValues = []  # we start with blanks for the values
    fieldValues = multenterbox(msg,title, fieldNames)
    print "Reply was:", fieldValues
if choice == 'p':
    multenterbox(msg='Fill in values for the fields.', title=' ', fields=(), values=())
    msg = "Enter your personal information"
    title = "Credit Card Application"
    fieldNames = ["Name","Street Address","City","State","ZipCode"]
    fieldValues = []  # we start with blanks for the values
    fieldValues = multenterbox(msg,title, fieldNames)
    print "Reply was:", fieldValues
if choice == 'e':
    multenterbox(msg='Fill in values for the fields.', title=' ', fields=(), values=())
    msg = "Enter your personal information"
    title = "Credit Card Application"
    fieldNames = ["Name","Street Address","City","State","ZipCode"]
    fieldValues = []  # we start with blanks for the values
    fieldValues = multenterbox(msg,title, fieldNames)
    print "Reply was:", fieldValues
if choice == 'randst':
    multenterbox(msg='Fill in values for the fields.', title=' ', fields=(), values=())
    msg = "Enter your personal information"
    title = "Credit Card Application"
    fieldNames = ["Name","Street Address","City","State","ZipCode"]
    fieldValues = []  # we start with blanks for the values
    fieldValues = multenterbox(msg,title, fieldNames)
    print "Reply was:", fieldValues
if choice == 'c':
    multenterbox(msg='Fill in values for the fields.', title=' ', fields=(), values=())
    msg = "Enter your personal information"
    title = "Credit Card Application"
    fieldNames = ["Name","Street Address","City","State","ZipCode"]
    fieldValues = []  # we start with blanks for the values
    fieldValues = multenterbox(msg,title, fieldNames)
    print "Reply was:", fieldValues
if choice == 'b':
    multenterbox(msg='Fill in values for the fields.', title=' ', fields=(), values=())
    msg = "Enter your personal information"
    title = "Credit Card Application"
    fieldNames = ["Name","Street Address","City","State","ZipCode"]
    fieldValues = []  # we start with blanks for the values
    fieldValues = multenterbox(msg,title, fieldNames)
    print "Reply was:", fieldValues
