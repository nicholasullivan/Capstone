import tkinter as tk
from tkinter import filedialog, Text
import os
from neighbor_joining2 import *

class App:

    def __init__(self, master):

        canvas=tk.Canvas(root, height=700, width=700, background='orange')
        canvas.pack()

        frame= tk.Frame(root, background='white')
        frame.place(relwidth=.9,relheight=.9,relx=.05,rely=.05)

        title=tk.Label(frame, text="DNA Sequence Similarity Calculator")

        start=tk.Button(frame,text="Start",padx=40, pady=20,command= lambda: Calculations.start())
        start.place(x= 275,y=375)

root = Tk()
app = App(root)

root.mainloop()