#! /usr/bin/env python

import sys
if sys.version_info > (2, 9):
  exit("must use python 2.7 or smaller due to lack of PIL module")

from Tkinter import *
from PIL.ImageTk import *

import os.path
import fileinput
import glob
import math
import string

# new development (2.4)  
#import gtk
import PIL.ImageDraw


######################
#  Global Variables  #
######################
TITLE   = "Ani2D, version 3.0"
aniLOGO = "logo.png"
aniNULL = "empty.png"

aniMESHini  = "../bin/mesh_initial.ps"
aniMESHfin  = "../bin/mesh_final.ps"
aniISOLINES = "../bin/iso_fin.ps"
aniVELOCITY = "../bin/streamlines_fin.ps"
aniPRESSURE = "../bin/pressure_fin.ps"
aniMATRIX   = "../bin/matrix_fin.ps"

pyMESHini = "../bin/mesh_initial.py"
pyMESHfin = "../bin/mesh_final.py"

meshINI = aniLOGO
meshFIN = aniLOGO


##################
#  Call AniDemo  #
##################

class AniDemo:
  def __init__(s, root):
    s.root = root

    s.varF1 = 0
    s.varG1 = 0
    s.varH1 = 0

    s.file  = aniLOGO
    s.file1 = aniNULL
    s.file2 = aniNULL

    #TOP MENU
    s.mBar = Frame(s.root, relief=RAISED, borderwidth=2)
    s.mBar.pack(fill=X)

    s.interpolation = s.interpolationMenu(s.mBar)
    s.packages      = s.packagesMenu(s.mBar)
    s.tutorial      = s.tutorialMenu(s.mBar)
    s.analysis      = s.analysisMenu(s.mBar)
    s.help          = s.helpMenu(s.mBar)

    #LEFT FRAME
    s.leftframe = Frame(s.root)
    s.leftframe.pack(side=LEFT)

    s.title0 = Label(s.leftframe, text="CLICK on MESH\n\n" \
                                       " or \n\n SELECT EXAMPLE\n\n\n", fg='dark blue')
    s.title0.pack()

    s.title1 = Label(s.leftframe, fg='dark blue')  # 1
    s.title1.pack()
    s.title1.config(text="Initial Mesh")

    s.img1 = Image.open(s.file1)
    s.img1.thumbnail((128, 128), Image.ANTIALIAS)
    s.photo1  = PhotoImage(s.img1)
    s.button1 = Button(s.leftframe,
                       background='gray', activebackground='dark red', borderwidth=3,
                       image=s.photo1, command=s.enlargeImage1)
    s.button1.pack(fill=X)

    s.title2 = Label(s.leftframe, fg='dark blue')  # 2
    s.title2.pack()
    s.title2.config(text="\nFinal Mesh")

    s.img2 = Image.open(s.file2)
    s.img2.thumbnail((128, 128), Image.ANTIALIAS)
    s.photo2  = PhotoImage(s.img2)
    s.button2 = Button(s.leftframe,
                       background='gray', activebackground='dark red', borderwidth=3,
                       image=s.photo2, command=s.enlargeImage2)
    s.button2.pack(fill=X)

    #RIGHT FRAME
    s.rightframe = Frame(s.root)
    s.rightframe.pack(side=RIGHT)

    s.title = Label(s.rightframe, fg='dark blue')
    s.title.pack()
    s.title.config(text="")

    s.rightcanvas = Canvas(s.rightframe, borderwidth=1)
    s.rightcanvas.pack()

    s.img = Image.open(s.file)
    s.photo = PhotoImage(s.img)
#   s.rightcanvas.create_image(0, 0, image=s.photo, anchor=Tkinter.NW)
    s.rightcanvas.create_image(0, 0, image=s.photo, anchor=NW)
    s.rightcanvas.config(width=s.img.size[0], height=s.img.size[1])
  
    # CALL INFORMATION WINDOW
    nw = Toplevel()
    nw.transient(root)
    nw.geometry("+200+120")
    nw.title("WARNING!")
    frame = Frame(nw, height=12, width=32)
    frame.pack()

    text = "MAKE SURE THAT THE PACKAGE HAS BEEN INSTALLED" 

    nw.label = Button(frame, bg='red', activebackground='white',
                      padx=5, pady=20, text=text, command=nw.destroy)
    nw.label.pack()


##################
#  Main program  #

  #######################################################
  #   MENUS
  #######################################################
  def interpolationMenu(s, mBar):
     """ Creates first menu list for package aniMBA 
     """
     s.ADAPT = IntVar()

     button = Menubutton(mBar, text='MESH ADAPTATION')
     button.pack(side=LEFT, padx='2m')

     # the primary pulldown
     button.menu = Menu(button)

     button.menu.add_radiobutton(label='Isotropic Metric',     
                                 value=1, variable=s.ADAPT, command=s.updateAdaptation)
     button.menu.add_separator()
     button.menu.add_radiobutton(label='Anisotropic Metric 1',       
                                 value=2, variable=s.ADAPT, command=s.updateAdaptation)
     button.menu.add_radiobutton(label='Stronger Anisotropic Metric 1', 
                                 value=3, variable=s.ADAPT, command=s.updateAdaptation)
     button.menu.add_separator()
     button.menu.add_radiobutton(label='Anisotropic Metric 2',       
                                 value=4, variable=s.ADAPT, command=s.updateAdaptation)
     button.menu.add_radiobutton(label='Stronger Anisotropic Metric 2',       
                                 value=5, variable=s.ADAPT, command=s.updateAdaptation)
     button.menu.add_separator()
     button.menu.add_radiobutton(label='Anisotropic Metric 3',       
                                 value=6, variable=s.ADAPT, command=s.updateAdaptation)
     button.menu.add_radiobutton(label='Stronger Anisotropic Metric 3',       
                                 value=7, variable=s.ADAPT, command=s.updateAdaptation)

     # set up a pointer from the file menubutton back to the file menu
     button['menu'] = button.menu

     return button


  def packagesMenu(s, mBar):
     s.PACKAGE = IntVar()

     button = Menubutton(mBar, text='PACKAGEs')
     button.pack(side=LEFT, padx='2m')

     button.menu = Menu(button)
     button.menu.add_radiobutton(label='Package Ani2D-AFT', 
                                 value=1, variable=s.PACKAGE, command=s.testPackages)
     button.menu.add_radiobutton(label='Package Ani2D-FEM', 
                                 value=2, variable=s.PACKAGE, command=s.testPackages)
     button.menu.add_radiobutton(label='Package Ani2D-MBA', 
                                 value=3, variable=s.PACKAGE, command=s.testPackages)
     button.menu.add_radiobutton(label='Package Ani2D-RCB',
                                 value=4, variable=s.PACKAGE, command=s.testPackages)
     button.menu.add_radiobutton(label='MultiPackage Anisotropic Adapatation',
                                 value=5, variable=s.PACKAGE, command=s.testPackages)
     button.menu.add_radiobutton(label='MultiPackage Convection Diffusion',
                                 value=6, variable=s.PACKAGE, command=s.testPackages)
     button.menu.add_radiobutton(label='MultiPackage Navier Stokes',
                                 value=7, variable=s.PACKAGE, command=s.testPackages)

     button['menu'] = button.menu

     return button


  def tutorialMenu(s, mBar):
     button = Menubutton(mBar, text='TUTORIALs')
     button.pack(side=LEFT, padx='2m')

     button.menu = Menu(button)
     button.menu.add_radiobutton(label='Interpolation problem', command=s.updateInterpolation)
     button.menu.add_command(label='Stokes problem', command=s.Stokes)

     button['menu'] = button.menu

     return button


  def analysisMenu(s, mBar):
     s.ANALYSIS = IntVar()

     button = Menubutton(mBar, text='ANALYSIS')
     button.pack(side=LEFT, padx='2m')

     button.menu = Menu(button)
     s.MENU = button.menu
     button.menu.add_radiobutton(label='Show final mesh', state=DISABLED,
                                 value=1, variable=s.ANALYSIS, command=s.updateAnalysis)

     button.menu.add_radiobutton(label='Show solution isolines', state=DISABLED,
                                 value=2, variable=s.ANALYSIS, command=s.updateAnalysis)

     button.menu.add_radiobutton(label='Show velocity streamlines', state=DISABLED,
                                 value=3, variable=s.ANALYSIS, command=s.updateAnalysis)

     button.menu.add_radiobutton(label='Pressure isolines', state=DISABLED,
                                 value=4, variable=s.ANALYSIS, command=s.updateAnalysis)

     button.menu.add_radiobutton(label='Matrix structure', state=DISABLED,
                                 value=5, variable=s.ANALYSIS, command=s.updateAnalysis)

     button['menu'] = button.menu
     return button


  def helpMenu(s, mBar):
     button = Menubutton(mBar, text='Help', underline=0)
     button.pack(side=RIGHT, padx='2m')

     button.menu = Menu(button)
     button.menu.add_command(label="About...", command=s.about)
     button.menu.add_separator()
     button.menu.add_command(label="Quit", command=root.quit)

     button['menu'] = button.menu
     return button


  #######################################################
  #  SUBROUTINES for IMAGES
  #######################################################
  def enlargeImage1(s):
     s.file = s.file1 
     s.title.config(text="Initial mesh")

     s.img   = Image.open(s.file)
     s.photo = PhotoImage(s.img)
#    s.rightcanvas.create_image(0, 0, image=s.photo, anchor=tkinter.NW)
     s.rightcanvas.create_image(0, 0, image=s.photo, anchor=NW)

  def enlargeImage2(s):
     s.file = s.file2
     s.title.config(text="Final quasi-uniform mesh")

     s.img   = Image.open(s.file)
     s.photo = PhotoImage(s.img)
#    s.rightcanvas.create_image(0, 0, image=s.photo, anchor=tkinter.NW)
     s.rightcanvas.create_image(0, 0, image=s.photo, anchor=NW)

  def updateImages(s):
     s.img1 = Image.open(s.file1)
     s.img1.thumbnail((128, 128), Image.ANTIALIAS)
     s.photo1  = PhotoImage(s.img1)
     s.button1.config(image=s.photo1)

     s.img2 = Image.open(s.file2)
     #s.img2 = s.drawMesh(pyMESHfin)
     s.img2.thumbnail((128, 128), Image.ANTIALIAS)
     s.photo2  = PhotoImage(s.img2)
     s.button2.config(image=s.photo2)

     s.img = Image.open(s.file)
     #s.img = s.drawMesh(pyMESHfin)
     s.photo = PhotoImage(s.img)
#    s.rightcanvas.create_image(0, 0, image=s.photo, anchor=tkinter.NW)
     s.rightcanvas.create_image(0, 0, image=s.photo, anchor=NW)

  def drawMesh(s, filename):
     img = Image.open("empty.png")
     draw = PIL.ImageDraw.Draw(img)

     f = open(filename)
     nt = int(f.readline())
     
     for i in range(1,nt):
        line = f.readline()
        objs = line.split()
        
        x1 = tuple([eval(objs[0]), eval(objs[1])])
        x2 = tuple([eval(objs[2]), eval(objs[3])])
        x3 = tuple([eval(objs[4]), eval(objs[5])])
      
        #draw.polygon([x1, x2, x3])
        draw.line(x1 + x2, fill=128)
        draw.line(x2 + x3, fill=128)
        draw.line(x3 + x1, fill=128)

     return img


  #######################################################
  #  SUBROUTINES: ADAPTATION
  #######################################################
  def updateAdaptation(s):
     i = s.ADAPT.get() 

     wait = Toplevel()
     wait.transient(root)
     wait.title("Wait")
     wait.grab_set()
     Label(wait, text="Adapting mesh, please wait...", 
                 padx=5, pady=5, 
                 fg="red", borderwidth=5).pack()

     s.root.config(cursor="watch")
     s.root.update_idletasks()

     #set up the tensor
     if i == 1:
       s.varF1 = 1 
       s.varG1 = 1 
       s.varH1 = 0
       s.title.config(text="Final quasi-optimal mesh for constant isotropic metric")
     elif i == 2:
       s.varF1 = 50 
       s.varG1 = 1 
       s.varH1 = 0
       s.title.config(text="Final quasi-optimal mesh for anisotropic metric [50 0; 0 1]")
     elif i == 3:
       s.varF1 = 200 
       s.varG1 = 1 
       s.varH1 = 0
       s.title.config(text="Final quasi-optimal mesh for anisotropic metric [200 0; 0 1]")
     elif i == 4:
       s.varF1 = 1 
       s.varG1 = 50 
       s.varH1 = 0
       s.title.config(text="Final quasi-optimal mesh for anisotropic metric [1 0; 0 50]")
     elif i == 5:
       s.varF1 = 1 
       s.varG1 = 200 
       s.varH1 = 0
       s.title.config(text="Final quasi-optimal mesh for anisotropic metric [1 0; 0 200]")
     elif i == 6:
       s.varF1 = 33.3799334 
       s.varG1 = 33.3799334
       s.varH1 =-32.0709164
       s.title.config(text="Final quasi-optimal mesh for anisotropic metric [33.4 -32.1; -32.1 33.4]")
     elif i == 7:
       s.varF1 = 131.556208
       s.varG1 = 131.556208
       s.varH1 =-130.247191
       s.title.config(text="Final quasi-optimal mesh for anisotropic metric [131.6 -130.2; -130.2 131.6]")

     #creare a file and solve the problem
     s.printAdaptation()
     os.system("cd ../bin; ./aniPY_mba.exe > /dev/null")

     s.file  = aniMESHfin
     s.file1 = aniMESHini
     s.file2 = aniMESHfin
     s.updateImages()

     wait.destroy()
     s.root.config(cursor="")

     s.MENU.entryconfigure(1, state="normal")
     s.MENU.entryconfigure(2, state="disabled")
     s.MENU.entryconfigure(3, state="disabled")
     s.MENU.entryconfigure(4, state="disabled")
     s.MENU.entryconfigure(5, state="disabled")


  def printAdaptation(s):
     file = open("../bin/aniPY_mba.txt", "w")
     file.write("%2i varF1\n" % s.varF1)
     file.write("%2i varG1\n" % s.varG1)
     file.write("%2i varH1\n" % s.varH1)
     file.close()


  #######################################################
  #  SUBROUTINES: PACKAGE AFT, MBA, FEM
  #######################################################
  def testPackages(s):
     i = s.PACKAGE.get() 

     wait = Toplevel()
     wait.transient(root)
     wait.title("Wait")
     wait.grab_set()

     Label(wait, text="\nTesting the selected Package or MultiPackage\n\n"
                      "MultiPackage may take couple of minutes\n"
                      "pay attention to the screen output...  \n\n"
                      "Initial and final mesh will be available for viewing\n",
                      padx=5, pady=5, 
                      fg="blue", borderwidth=5).pack()

     s.root.config(cursor="watch")
     s.root.update_idletasks()

     if i == 1:
        os.system("cd ../bin; ./aniAFT_front1.exe")
     elif i == 2:
        os.system("cd ../bin; ./aniFEM_template.exe")
     elif i == 3:
        os.system("cd ../bin; ./aniMBA_analytic.exe")
     elif i == 4:
        os.system("cd ../bin; ./aniRCB.exe")
     elif i == 5:
        os.system("cd ../bin; ./multi_AnisoAdap2.exe")
     elif i == 6:
        os.system("cd ../bin; ./multi_ConDifAdap.exe")
     elif i == 7:
        os.system("cd ../bin; ./multi_StokesNavier.exe")

     wait.bell()
     s.root.config(cursor="")

     meshFIN = aniMESHfin
     s.file1 = aniMESHini
     s.file2 = aniMESHfin
     s.updateImages()
     wait.destroy()

     for k in range(1, 6):
         s.MENU.entryconfigure(k, state="disabled")

     if i == 1:
        s.MENU.entryconfigure(1, state="normal")
     elif i == 2:
        s.MENU.entryconfigure(5, state="normal")
     elif i == 3:
        s.MENU.entryconfigure(1, state="normal")
     elif i == 4:
        s.MENU.entryconfigure(1, state="normal")
     elif i > 4:
        s.MENU.entryconfigure(1, state="normal")
        s.MENU.entryconfigure(2, state="normal")
        s.MENU.entryconfigure(5, state="normal")

     if i == 7: 
        s.MENU.entryconfigure(3, state="normal")


  #######################################################
  #  SUBROUTINES: ITERPOLATION
  #######################################################
  def updateInterpolation(s):
     wait = Toplevel()
     wait.transient(root)
     wait.title("Wait")
     wait.grab_set()

     Label(wait, text="Solving problem, please wait...", fg="red", borderwidth=5).pack()

     s.root.config(cursor="watch")
     s.root.update_idletasks()

     os.system("cd ../bin; ./multi_Interpolation.exe > /dev/null")
     s.title.config(text="Final quasi-optimal mesh minimizing the maximum norm\n" \
                         "of the interpolation error for piecewise linear approximations")

     s.file  = aniMESHfin
     s.file1 = aniMESHini
     s.file2 = aniMESHfin
     s.updateImages()

     wait.bell()
     wait.destroy()
     s.root.config(cursor="")

     s.MENU.entryconfigure(1, state="normal")
     s.MENU.entryconfigure(2, state="normal")
     s.MENU.entryconfigure(3, state="disabled")
     s.MENU.entryconfigure(4, state="disabled")
     s.MENU.entryconfigure(5, state="disabled")


  #######################################################
  #  SUBROUTINES: STOKES
  #######################################################
  def Stokes(s):
     nw = Toplevel()
     nw.transient(root)
     nw.title("Adaptive solution of the Stokes problem")

     # define boundary conditions (Left and Right)
     domain = Label(nw, text="Boundary Conditions", fg='dark blue')
     domain.grid(row=3, column=0, sticky=W, columnspan=3, padx=10, pady=0)

     s.BC_LEFT = StringVar()
     domain = Label(nw, text="BC LEFT:")
     domain.grid(row=4, column=0, sticky=W, padx=10, pady=0)
     button = OptionMenu(nw, s.BC_LEFT, "Dirichlet", "Neumann")
     button.grid(row=4, column=1, sticky=W, padx=10, pady=0)

     s.BC_RIGHT = StringVar()
     domain = Label(nw, text="BC RIGHT:")
     domain.grid(row=4, column=2, sticky=W, padx=10, pady=0)
     button = OptionMenu(nw, s.BC_RIGHT, "Dirichlet", "Neumann")
     button.grid(row=4, column=3, sticky=W, padx=10, pady=0)

     # define boundary conditions (Top and Bottom)
     s.BC_TOP = StringVar()
     domain = Label(nw, text="BC TOP:")
     domain.grid(row=5, column=0, sticky=W, padx=10, pady=0)
     button = OptionMenu(nw, s.BC_TOP, "No Slip", "Neumann")
     button.grid(row=5, column=1, sticky=W, padx=10, pady=0)

     s.BC_BOTTOM = StringVar()
     domain = Label(nw, text="BC BOTTOM:")
     domain.grid(row=5, column=2, sticky=W, padx=10, pady=0)
     button = OptionMenu(nw, s.BC_BOTTOM, "No Slip", "Neumann")
     button.grid(row=5, column=3, sticky=W, padx=10, pady=0)

     # define source term
     domain = Label(nw, text="\nSource Term", fg='dark blue')
     domain.grid(row=6, column=0, sticky=W, columnspan=2, padx=10, pady=0)

     s.SOURCEX = StringVar()
     domain = Label(nw, text="Function:")
     domain.grid(row=7, column=0, sticky=W, padx=10, pady=0)
     button = OptionMenu(nw, s.SOURCEX, "Q_x=0", "Q_x=x")
     button.grid(row=7, column=1, sticky=W, padx=10, pady=0)

     # define finite elements
     domain = Label(nw, text="\nSelect Finite Elements, unstable pair is allowed for FUN", fg='dark red')
     domain.grid(row=8, column=0, sticky=W, columnspan=3, padx=10, pady=0)

     s.FEM_U = StringVar()
     domain = Label(nw, text="FEM VELOCITY:")
     domain.grid(row=9, column=0, sticky=W, padx=10, pady=0)
     button = OptionMenu(nw, s.FEM_U, "FEM_P2vector", "FEM_MINI", "FEM_P2reduced")
     button.grid(row=9, column=1, sticky=W, padx=10, pady=0)

     s.FEM_P = StringVar()
     domain = Label(nw, text="FEM PRESSURE:")
     domain.grid(row=10, column=0, sticky=W, padx=10, pady=0)
     button = OptionMenu(nw, s.FEM_P, "FEM_P1", "FEM_P0")
     button.grid(row=10, column=1, sticky=W, padx=10, pady=0)

     # define finite elements
     domain = Label(nw, text="\nAdaptation Parameters", fg='dark blue')
     domain.grid(row=11, column=0, sticky=W, columnspan=2, padx=10, pady=0)

     s.PHYSICS = StringVar()
     domain = Label(nw, text="STRATEGY")
     domain.grid(row=12, column=0, sticky=W, padx=10, pady=0)
     button = OptionMenu(nw, s.PHYSICS, "Adaptation to U_x", "Adaptation to U_y")
     button.grid(row=12, column=1, columnspan=3, sticky=W, padx=10, pady=0)

     s.QUALITY = DoubleVar()
     label = Label(nw, text="MESH QUALITY")
     label.grid(row=13, column=0, sticky=W, padx=10, pady=0)
     button = OptionMenu(nw, s.QUALITY, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)
     button.grid(row=13, column=1, sticky=W, padx=10, pady=0)

     s.NESTAR = IntVar()
     domain = Label(nw, text="FINAL NUMBER of TRIANGLES")
     domain.grid(row=14, column=0, columnspan=2, sticky=W, padx=10, pady=0)
     button = OptionMenu(nw, s.NESTAR, 5000, 10000, 20000)
     button.grid(row=14, column=2, sticky=W, padx=10, pady=0)

     s.NLOOPS = IntVar()
     label = Label(nw, text="TOTAL NUMBER of ADAPTATION LOOPS")
     label.grid(row=15, column=0, columnspan=2, sticky=W, padx=10, pady=0)
     button = OptionMenu(nw, s.NLOOPS, 1, 2, 3, 4, 5)
     button.grid(row=15, column=2, sticky=W, padx=10, pady=0)

     # initial values
     s.BC_LEFT.set("Dirichlet")
     s.BC_RIGHT.set("Neumann")
     s.BC_TOP.set("No Slip")
     s.BC_BOTTOM.set("No Slip")

     s.SOURCEX.set("Q_x=0")

     s.FEM_U.set("FEM_P2vector")
     s.FEM_P.set("FEM_P1")

     s.PHYSICS.set("Adaptation to U_x")
     s.QUALITY.set(0.7)
     s.NESTAR.set(5000)
     s.NLOOPS.set(2)

     # executable
     label = Label(nw, text="\n-------------------------------------------------" \
                  "-------------------------------------------------------------", fg='dark blue')
     label.grid(row=16, column=0, columnspan=4, sticky=W, padx=10, pady=0)
     label = Button(nw, activebackground='white', text="Save and Execute", command=s.updateStokes)
     label.grid(row=17, column=1, columnspan=2, padx=10, sticky=S, pady=0)


  def updateStokes(s):
     wait = Toplevel()
     wait.transient(root)
     wait.title("Wait")
     wait.grab_set()
     wait.focus_set()

     Label(wait, text="Solving the Stokes problem, please wait...", fg="red", borderwidth=5).pack()

     s.root.config(cursor="watch")
     s.root.update_idletasks()

     #creare a file and solve the problem
     s.printStokes()
     os.system("cd ../bin; ./aniPY_stokes.exe > /dev/null")

     s.file1 = aniMESHini
     s.file2 = aniMESHfin
     s.updateImages()

     wait.bell()
     wait.destroy()
     s.root.config(cursor="")

     s.MENU.entryconfigure(1, state="normal")
     s.MENU.entryconfigure(2, state="disabled")
     s.MENU.entryconfigure(3, state="normal")
     s.MENU.entryconfigure(4, state="normal")
     s.MENU.entryconfigure(5, state="normal")


  def printStokes(s):
     phys = s.PHYSICS.get()
     if phys == "Adaptation to U_x":
        phys_id = 0
     elif phys == "Adaptation to U_y":
        phys_id = 1

     bc = s.BC_TOP.get()
     if bc == "No Slip":
        bc_top = 1
     elif bc == "Neumann":
        bc_top = 0

     bc = s.BC_BOTTOM.get()
     sourceText  = "PY_BC_BOTTOM"
     if bc == "No Slip":
        bc_bottom = 1
     elif bc == "Neumann":
        bc_bottom = 0

     bc = s.BC_LEFT.get()
     if bc == "Dirichlet":
        bc_left = 1
     elif bc == "Neumann":
        bc_left = 0

     bc = s.BC_RIGHT.get()
     if bc == "Dirichlet":
        bc_right = 1
     elif bc == "Neumann":
        bc_right = 0

     src = s.SOURCEX.get()
     if src == "Q_x=0":
        src_id = 0
     elif src == "Q_x=x":
        src_id = 1

     fem = s.FEM_U.get()
     if fem == "FEM_P2vector":
        fem_u = 12
     elif fem == "FEM_P2reduced":
        fem_u = 13
     elif fem == "FEM_MINI":
        fem_u = 14

     fem = s.FEM_P.get()
     if fem == "FEM_P1":
        fem_p = 2
     elif fem == "FEM_P0":
        fem_p = 1

     file = open("../bin/aniPY_stokes.txt", "w")
     file.write("%7d   NLOOPS\n" % s.NLOOPS.get())
     file.write("%7d   nEstar\n" % s.NESTAR.get())
     file.write("%6.2f Quality\n" % s.QUALITY.get())
     file.write("%7d   Physics\n" % phys_id)
     file.write("%7d   BC_TOP\n" % bc_top)
     file.write("%7d   BC_BOTTOM\n" % bc_bottom)
     file.write("%7d   BC_LEFT\n" % bc_left)
     file.write("%7d   BC_RIGHT\n" % bc_right)
     file.write("%7d   SOURCE_X\n" % src_id)
     file.write("%7d   FEM_U\n" % fem_u)
     file.write("%7d   FEM_P\n" % fem_p)
     file.close()


  #######################################################
  #  SUBROUTINES: ANALYSIS
  #######################################################
  def updateAnalysis(s):
     i = s.ANALYSIS.get() 

     if i == 1:
       s.title.config(text="Final quasi-optimal mesh")
       s.file  = aniMESHfin
     elif i == 2:
       s.title.config(text="Solution isolines")
       s.file  = aniISOLINES
     elif i == 3:
       s.title.config(text="Velocity streamlines")
       s.file  = aniVELOCITY
     elif i == 4:
       s.title.config(text="Pressure isolines")
       s.file  = aniPRESSURE
     elif i == 5:
       s.title.config(text="Matrix sparcity structure")
       s.file  = aniMATRIX

     s.updateImages()


  #######################################################
  #  SUBROUTINES: HELP
  #######################################################
  def about(self):
     nw = Toplevel()
     nw.transient(root)
     nw.title("About")
     frame = Frame(nw, height=12, width=32)
     frame.pack()

     text = "MESH ADAPTATION. Select a metric from\n"   \
            "the first menu and find out which mesh\n"  \
            "is quasi-uniform in this metric\n\n"       \
            "PACKAGEs. Capabilities of individual  \n"  \
            "packages as well a few multi-packages.\n"   \
            "Package Ani2D has been installed first\n\n" \
            "TUTORIALs. Examples of solving complex\n"  \
            "problems: optimal mesh for minimizing\n"   \
            "the interpolation error, finite element\n" \
            "solution of the Stokes problem, etc.\n\n"  \
            "ANALYSIS. View details of the complex\n"   \
            "examples: meshes, solution isolines,\n"    \
            "streamlines, etc.\n\n\n"                   \
            "DEVELOPER: Konstantin Lipnikov\n"          \
            "\a lipnikov@gmail.com"

     nw.label = Button(frame, activebackground='white',
                       padx=5, pady=30, text=text, command=nw.destroy)
     nw.label.pack()



##################
#  Main program  #
##################

root = Tk()
root.wm_geometry("660x580+20+20")
root.title(TITLE)

demo = AniDemo(root)

root.mainloop()


