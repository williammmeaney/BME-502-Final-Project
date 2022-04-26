#Protein Calc
#ug/ul=Fqn*50
#Protein=960/ug
#CL=204-Protein
using CSV# as CV
using DataFrames
using Gtk
#Import
#file_name=open_dialog("Pick a file");
#file_name
file_contents=CSV.read("C:\\Users\\willi\\Downloads\\DataForFinal.csv", DataFrame);
file_contents_matrix=Matrix(file_contents);
#Extract Arrays in Float64
Absorbance=file_contents_matrix[:,2];
Abs = Array{Float64}(Absorbance)
B_actin=file_contents_matrix[:,3];
B_actin_Fl = Array{Float64}(B_actin)
Pink_1=file_contents_matrix[:,4];
Pink_1_Fl=Array{Float64}(Pink_1)
#Protein Calculations
Frequency=(Abs.-0.0915)./ 0.5759;
Concentration=Frequency.*50;
Protien=960 ./Concentration;
CL=204 .- Protien;
#Western Blot
max_b_actin=maximum(B_actin_Fl);
Corrected_B_actin=max_b_actin./ B_actin_Fl;
Pink_1_Corr=Corrected_B_actin ./ Pink_1_Fl;

#corrected B actin - remove above 3 for WB only

