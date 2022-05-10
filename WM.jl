#Protein Calc
#ug/ul=Fqn*50
#Protein=960/ug
#CL=204-Protein
using CSV# as CV
using DataFrames
using Gtk
using Statistics
#Import
#file_name=open_dialog("Pick a file");

#file_name
file_contents=CSV.read("C:\\Users\\willi\\Downloads\\DataForFinal.csv", DataFrame);
file_contents_matrix=Matrix(file_contents);
#Extract Arrays in Float64
Groups=file_contents_matrix[:,1];
list_of_unique_groups=unique(Groups);
Absorbance=file_contents_matrix[:,3];
Abs = Array{Float64}(Absorbance)
B_actin=file_contents_matrix[:,4];
B_actin_Fl = Array{Float64}(B_actin)
Pink_1=file_contents_matrix[:,5];
Pink_1_Fl=Array{Float64}(Pink_1)

#groups 

#Grouping
#=#number in each group is stated upfront in fucntion ahead of time
#validation 
if length(unique_group_inputs)==length(list_of_unique_groups)
    number_groups_match=true
else
    number_groups_match=false
end
=#


#create array of begining and end rows for each group
#number_in_group - holds number of samples in each Group
# input_group_names - holds names of groups
group_values=zeros(length(unique_group_inputs),2)
group_values(1,1)=1;
group_values(2,2)=number_in_group(1);
counter=number_in_group(1);
for i=2:length(input_group_names)
    counter=counter+1;
    group_values(i,1)=counter;
    group_values(i,2)=counter+number_in_group(i+1)-1;
end


#Determine unique categry names and number of elements in each category
unique_categories=unique(input_group_names);
number_unique=length(unique_categories);
count_of_groups=ones(number_unique);
for i=1:number_unique
    count_of_groups(i)=count(==(unique_categories(i)), input_group_names)
end

#manual input
#=
print("What is the total number of groups you have (number experimental and control)\n\n");
total=readline();
input_group_names=[];
number_in_group=[];
print("You will be prompted for several group names and number of samples in each group one by one. Please go in order of the csv file you inputted. \n\n")
for i = 1:total
   print("Please state a group name.\n\n") 
   group_name_input=readline();
   push!(input_group_names,group_name_input)
   print("Please state the number of samples in each group.\n\n") 
   number_group_input=readline();
   push!(number_in_group,nunber_group_input)
end
=#



#Protein Calculations
Frequency=(Abs.-0.0915)./ 0.5759;
Concentration=Frequency.*50;
Protien=960 ./Concentration;
CL=204 .- Protien;
#Western Blot
max_b_actin=maximum(B_actin_Fl);


#corrected B actin - remove above 3 for WB only
B_actin_Fl_below_3=[];
Pink_1_Corr_below_3=[];
Groups_below_3=[];
show_outliers=false;
if show_outliers=true #add to function
    Corrected_B_actin_w_outliers=max_b_actin./ B_actin_Fl;
    Pink_1_Corr_w_outliers=Corrected_B_actin ./ Pink_1_Fl;
else 
    for i=1:length(B_actin_Fl)
        if B_actin_Fl(i)<3
            push!(B_actin_Fl_below_3,B_actin_Fl(i));
            push!(Pink_1_Fl_below_3,Pink_1_Fl(i));
            push!(Groups_below_3,Groups(i));
        end
    end
    Corrected_B_actin=max_b_actin./ B_actin_Fl_below_3;
    Pink_1_Corr=Corrected_B_actin_below_3 ./ Pink_1_Fl_below_3;
end





#gausean 
#mean and std --> percentile of (z=3) --> how many data points outside of percentile --> recalculate without outlier and 
#recaclualte mean and STD and then see how likley it would be outside

#measure how much outside by using gausean didstribution

g_mean=[];
g_sd=[]
for i=1:number_unique
    push!(g_mean,average(B_actin_Fl[group_values[i,1],group_values[i,2]]))
    push!(g_sd,stdm(g_mean,B_actin_Fl[group_values[i,1],group_values[i,2]]))

end