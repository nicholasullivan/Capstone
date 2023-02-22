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
        self.title=tk.Label(self.frame, text="Protein Sequence\nDistance Calculator",font= ("lucida 20 bold italic", 16))
        self.title.place(x=230,y=30)

        #add image
        self.pic=Image.open('capstone_image.jpg')
        self.pic=self.pic.resize((200, 200))
        self.pic=ImageTk.PhotoImage(self.pic)
        self.pic_win = tk.Label(root, image = self.pic)
        self.pic_win.place(x=33,y=33)

        #select file button
        self.file=tk.Button(self.frame,text="File Select",font= ("lucida 20 bold italic", 10),padx=20, pady=10, command = lambda:[Calculations.file_select(), enable_submit()])
        self.file.place(x= 270,y=180)

        '''#error popup window if wrong file is inputted
        def error_popup():
            self.top = tk.Toplevel(self.frame)
            self.top.geometry()
            self.top.geometry("300x200")
            self.top.title("File Type Error")
            tk.Label(self.top, text= 
                "Incompatible file type selected. Please choose a Phylip(.phy) or Fasta(.fa) msa file.",
                font=('lucida 12')).place(x=100,y=100)

        #file path printed
        def print_path():
            a = "File:" + str(Calculations.filename)
            self.path = tk.Label(self.frame, text=a)
            self.path.place(x=75, y=275)'''

        #quit button
        self.quit=tk.Button(self.frame, text="Quit App", font= ("lucida 20 bold italic", 10), padx=10, pady=10, command=root.destroy)
        self.quit.place(x= 540,y=18)

        #drop-down menu for score matrix selection
        matrices = [
            'Auto-assign BLOSUM based on identity','Auto-assign PAM based on identity',
            'BLOSUM30','BLOSUM40','BLOSUM50','BLOSUM62','BLOSUM70','BLOSUM80','BLOSUM90','BLOSUM100',
            'PAM10','PAM100','PAM200','PAM300','PAM400','PAM500'
        ]

        self.var = tk.StringVar()
        self.var.set('Choose a scoring matrix:') #default
        self.menu = tk.OptionMenu(self.frame, self.var, *matrices)
        self.menu.place(x = 230, y = 350)

        #submit button
        self.submit = tk.Button(self.frame, text="Submit", font= ("lucida 20 bold italic", 8), padx=8, pady=5, command=lambda:[Calculations.matrix_selection(self.var.get()),Calculations.show_tree('./nj.tree')], state = 'disabled')
        self.submit.place(x= 300,y=400)

        
        
        #function to enable button when file is selected
        def enable_submit():
            self.submit['state'] = 'normal'

        
root = Tk()
app = App(root)

root.mainloop()