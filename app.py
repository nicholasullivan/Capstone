import tkinter as tk
from tkinter import filedialog, Text
import os
from neighbor_joining2 import *

class App:

    def __init__( self, master):

        self.canvas=tk.Canvas(root, height=700, width=700, background='orange')
        self.canvas.pack()

        self.frame= tk.Frame(root, background='white')
        self.frame.place(relwidth=.9,relheight=.9,relx=.05,rely=.05)

        self.title=tk.Label(self.frame, text="DNA Sequence Similarity Calculator")

        self.start=tk.Button(self.frame,text="Start",font= ("lucida 20 bold italic", 14),padx=40, pady=20,command= lambda: Calculations.start())
        self.start.place(x= 270,y=375)

        self.quit=tk.Button(self.frame,text="Quit App",font= ("lucida 20 bold italic", 10),padx=10, pady=10,command= root.destroy)
        self.quit.place(x= 540,y=18)

root = Tk()
app = App(root)

root.mainloop()