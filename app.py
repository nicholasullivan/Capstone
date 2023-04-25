import tkinter as tk
from tkinter import Tk
from PIL import Image, ImageTk
from neighbor_joining2 import Calculations

class Interface:

    def __init__( self, master):
        #create a canvas
        self.canvas=tk.Canvas(root, height=650, width=700, background='#003DA5')
        self.canvas.pack()
        #create a frame
        self.frame= tk.Frame(root, background='white')
        self.frame.place(relwidth=.9,relheight=.85,relx=.05,rely=.03)

        #add title to frame
        self.title=tk.Label(self.frame, text="ProDisC",fg='#003DA5',bg='white',font=("Times New Roman", 30, "bold"))
        self.title.place(x=80,y=35)

        #add images
        self.tree = Image.open('assets/imgs/free_tree.png')
        width, height = self.tree.size
        self.tree = self.tree.resize((round(width*0.17), round(height*0.17)))
        self.tree = ImageTk.PhotoImage(self.tree)
        self.tree_win = tk.Label(root, image = self.tree)
        self.tree_win.place(x=350,y=40)

        self.logo = Image.open('assets/imgs/logo.png')
        width, height = self.logo.size
        self.logo = self.logo.resize((round(width*.12), round(height*.12)))
        self.logo = ImageTk.PhotoImage(self.logo)
        self.logo_win = tk.Label(root, image = self.logo)
        self.logo_win.place(x = round(700-width*.12)/2, y=585)

        #select file button
        file_lab = tk.Label(self.frame, text = 'Choose a Fasta or Phylip MSA file:', font=('lucida',10), background='white')
        file_lab.place(x=30, y=110)
        self.btn_text = tk.StringVar()
        self.btn_text.set("Select File")
        self.file_btn = tk.Button(self.frame, textvariable=self.btn_text, font=('lucida',10), padx=9, pady=5, command=lambda:[Calculations.file_select(), self.btn_text.set("File Selected")])
        
        self.file_btn.place(x=30, y=140)

        #drop-down menu for score matrix selection
        matrices = [
            'Auto-assign BLOSUM based on average identity score', 'Auto-assign BLOSUM based on pairwise identity score',
            'BLOSUM30','BLOSUM40','BLOSUM50','BLOSUM62','BLOSUM70','BLOSUM80','BLOSUM90',
            'PAM10','PAM100','PAM200','PAM300','PAM400','PAM500'
        ]
        
        mat_lab = tk.Label(self.frame, text = 'Choose a Scoring Matrix Option:', font=('lucida',10), background='white')
        mat_lab.place(x=30, y=180)
        self.mat = tk.StringVar()
        self.mat.set('')
        self.mat.menu = tk.OptionMenu(self.frame, self.mat, *matrices)
        self.mat.menu.place(x=30, y=210)

        #drop-down menu for gap penalty selection
        pens = ['-1','-2','-3','-4']
        gap_lab = tk.Label(self.frame, text = 'Choose a Gap Penalty:', font=('lucida',10), background='white')
        gap_lab.place(x=30, y=260)
        self.gap = tk.StringVar()
        self.gap.set('')
        self.gap.menu = tk.OptionMenu(self.frame, self.gap, *pens)
        self.gap.menu.place(x=30, y=290)

        #bootstrapping number of repetitions n entry box
        boot_lab = tk.Label(self.frame, text = 'Enter N number of bootstrap repetitions:', font=('lucida',10), background='white')
        boot_lab.place(x=30, y=330)
        self.n = tk.IntVar()
        n_entry = tk.Entry(self.frame, textvariable = self.n, font=('lucida',10))
        n_entry.place(x=30, y=360)

        #function to enable Run button when prerequisites are met
        self.file_warn = tk.Label(self.frame, text='No file was selected', font=('lucida',10), fg='white', bg='white')
        self.mat_warn = tk.Label(self.frame, text='No scoring matrix option was selected', font=('lucida',10), fg='white', bg='white')
        self.gap_warn = tk.Label(self.frame, text='No gap penalty was selected', font=('lucida',10), fg='white', bg='white')
        self.n_warn = tk.Label(self.frame, text='Number of bootstrap repititions must be greater than zero', font=('lucida',10), fg='white', bg='white')
        self.file_warn.place(x=30,y=450)
        self.mat_warn.place(x=30,y=470)
        self.gap_warn.place(x=30,y=490)
        self.n_warn.place(x=30,y=510)
        
        def check_readiness():
            self.file_warn.config(fg='white')
            self.mat_warn.config(fg='white')
            self.gap_warn.config(fg='white')
            self.n_warn.config(fg='white')
            if (self.btn_text.get()=='File Selected') & (self.n.get() > 0) & (self.mat.get() != ''):
                self.run['state'] = 'normal'
            if self.btn_text.get()=='Select File':
                self.file_warn.config(fg='red')
            if self.mat.get()=='':
                self.mat_warn.config(fg='red')
            if self.gap.get()=='':
                self.gap_warn.config(fg='red')
            if self.n.get()<=0:
                self.n_warn.config(fg='red')

        def popup(tf):
            while tf==True:
                if Calculations.complete()==True:
                    tk.messagebox.showinfo(title="Calculation Complete!", message="Consensus Tree and Summary files created\n in the same file location as MSA")
                break

        #submit buttons
        self.check = tk.Button(self.frame, text='Check', font= ("lucida", 16), padx=6, pady=2, command=lambda:check_readiness())
        self.run = tk.Button(self.frame, text="Run", font= ("lucida", 16), padx=10, pady=2,
                                command=lambda:[Calculations.matrix_selection(self.mat.get(), self.gap.get(),self.n.get()),popup(True)], state = 'disabled')
        self.check.place(x=40, y=390)
        self.run.place(x=165, y=390)

        #create and update progress label
        self.progress = tk.StringVar()
        self.prg_lab = tk.Label(self.frame, textvariable=self.progress, font=('lucida',10), fg='green', bg='white')
        
root = Tk()
root.title("Protein Sequence Distance Calculator")
app = Interface(root)

root.protocol("WM_DELETE_WINDOW", root.quit)
root.mainloop()