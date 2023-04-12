import tkinter as tk
from tkinter import filedialog, Text
import os
from neighbor_joining2 import *
from PIL import Image, ImageTk

class App:

    def __init__( self, master):
        #create a canvas
        self.canvas=tk.Canvas(root, height=700, width=700, background='orange')
        self.canvas.pack()
        #create a frame
        self.frame= tk.Frame(root, background='white')
        self.frame.place(relwidth=.9,relheight=.9,relx=.05,rely=.05)

        #add title to frame
        self.title=tk.Label(self.frame, text="Protein Sequence\nDistance Calculator",fg='orange',bg='white',font= ("lucida 20 bold italic", 18))
        self.title.place(x=210,y=30)

        #add image
        self.pic=Image.open('assets/imgs/capstone_image.jpg')
        self.pic=self.pic.resize((300, 300))
        self.pic=ImageTk.PhotoImage(self.pic)
        self.pic_win = tk.Label(root, image = self.pic)
        self.pic_win.place(x=200,y=170)

        #select file button
        self.file=tk.Button(self.frame, text = "MSA File Select", font= ("lucida 20 bold italic", 9), padx=9, pady=5, command=lambda:[Calculations.file_select(), enable_submit()])
        self.file.place(x= 170,y=460)

        #quit button
        self.quit=tk.Button(self.frame, text="Quit App", font= ("lucida 20 bold italic", 10),padx=10, pady=10, command=root.quit)
        self.quit.place(x= 540,y=18)

        #drop-down menu for score matrix selection
        matrices = [
            'Auto-assign BLOSUM based on average identity score', 'Auto-assign BLOSUM based on pairwise identity score'
            'BLOSUM30','BLOSUM40','BLOSUM50','BLOSUM62','BLOSUM70','BLOSUM80','BLOSUM90',
            'PAM10','PAM100','PAM200','PAM300','PAM400','PAM500'
        ]

        self.var = tk.StringVar()
        self.var.set('Choose a scoring matrix:') #default
        self.menu = tk.OptionMenu(self.frame, self.var, *matrices)
        self.menu.place(x=290, y=460)

        #bootstrapping number of repetitions n entry box
        n = tk.IntVar()
        n_lab = tk.Label(self.frame, text = 'N Bootstrap Reps', font=('lucida',10))
        n_entry = tk.Entry(self.frame, textvariable = n, font=('lucida',10))
        n.set(1)
        n_lab.place(x=180, y=510)
        n_entry.place(x=300, y=510)

        #submit button
        self.submit = tk.Button(self.frame, text="Submit", font= ("lucida 20 bold italic", 8), padx=8, pady=5,
                                command=lambda:[Calculations.matrix_selection(self.var.get(), n.get())], state = 'disabled')
        self.submit.place(x=290, y=560)

        #function to enable button when file is selected
        def enable_submit():
            self.submit['state'] = 'normal'
        
root = Tk()
root.title("Protein Sequence Distance Calculator")
app = App(root)

root.mainloop()