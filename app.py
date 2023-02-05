import tkinter as tk
from tkinter import filedialog, Text
import os
from neighbor_joining2 import *
from PIL import Image, ImageTk

class App:

    def __init__( self, master):

        self.canvas=tk.Canvas(root, height=700, width=700, background='orange')
        self.canvas.pack()

        self.frame= tk.Frame(root, background='white')
        self.frame.place(relwidth=.9,relheight=.9,relx=.05,rely=.05)

        """def open_new_win():
            top=tk.Toplevel(self.canvas)
            self.canvas2=tk.Canvas(root, height=180, width=100, bg="#aaaffe")
            self.canvas2.pack()
            tk.Label(self.canvas2, text="Tree File", font='Helvetica 18 bold').pack()"""

        self.title=tk.Label(self.frame, text="Protein Sequence\nDistance Calculator",font= ("lucida 20 bold italic", 16))
        self.title.place(x=230,y=30)

        self.pic=Image.open('capstone_image.jpg')
        self.pic=self.pic.resize((200, 200))
        self.pic=ImageTk.PhotoImage(self.pic)
        self.pic_win = tk.Label(root, image = self.pic)
        self.pic_win.place(x=33,y=33)


        self.file=tk.Button(self.frame,text="File Select",font= ("lucida 20 bold italic", 10),padx=20, pady=10,command= lambda: Calculations.file_select())
        self.file.place(x= 270,y=180)

        self.quit=tk.Button(self.frame,text="Quit App",font= ("lucida 20 bold italic", 10),padx=10, pady=10,command= root.destroy)
        self.quit.place(x= 540,y=18)

        self.matrix=tk.Label(self.frame, text='Scoring Matrix Options:\n'
		'1- BLOSUM\n'
		'2- PAM\n',padx=30, pady=15)
        self.matrix.place(x= 230,y=250)

        self.selection=tk.StringVar()
        self.mat=tk.Entry(self.frame, width=15,textvariable=self.selection)
        self.mat.place(x= 250,y=350)

        #self.file=tk.Button(self.frame,text="Get Tree File",font= ("lucida 20 bold italic", 10),padx=10, pady=10,command= self.canvas.text = canvas.create_text(20, 30, text=Calculations.get_file()))
        #self.file.place(x= 500,y=118)
        """def submission(text):
            Calculations.matrix_selection(text)"""
        self.submit=tk.Button(self.frame,text="Submit",font= ("lucida 20 bold italic", 8),padx=8, pady=5,command= lambda: Calculations.matrix_selection(self.mat.get()))
        self.submit.place(x= 350,y=350)

        
root = Tk()
app = App(root)

root.mainloop()