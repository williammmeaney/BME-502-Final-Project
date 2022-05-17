using CSV
using DataFrames
using Gtk
using Statistics
using Turing
using StatsPlots
using FillArrays
import SpecialFunctions as SF
import Distributions as Di
using HypothesisTests
import Pingouin as Pg

#=
The following cells before the csv import cells are attempts at metaprogramming - to dynamically create variables from the contents of a string
=#

#=function string_as_varname(s::String,v::Any)
         s=Symbol(s)
         @eval (($s) = ($v))
end
Symbol("sampl
@eval ($:sample
test_string="abc123"
:test_string
@eval($:string(abc123)=8)
=#

file_contents_vert=CSV.read("C:\\Users\\willi\\OneDrive - Stony Brook University\\Courses\\Undergraduate\\Senior Year\\BME 502\\GitHub\\BME-502-Final-Project\\DataForFinal_vert.csv", DataFrame);
file_contents_hor=CSV.read("C:\\Users\\willi\\OneDrive - Stony Brook University\\Courses\\Undergraduate\\Senior Year\\BME 502\\GitHub\\BME-502-Final-Project\\DataForFinal_hor.csv", DataFrame);
file_contents=copy(file_contents_vert)
file_contents_no_outliers=copy(file_contents_vert)


#protein calc
Frequency=(file_contents[:,:Absorbance].-0.0915)./ 0.5759;
Concentration=Frequency.*50;
Protien=960 ./Concentration;
CL=204 .- Protien;
#Western Blot
max_b_actin=maximum(file_contents[:,:B_actin]);

normalized_file_contents=copy(file_contents)
normalized_file_contents[:,:B_actin]=normalized_file_contents[:,:B_actin]/max_b_actin
normalized_file_contents[:,:Pink_1]=normalized_file_contents[:,:B_actin].*normalized_file_contents[:,:Pink_1]
grouped_trials=groupby(normalized_file_contents,:SampleType)

Pink_control=grouped_trials[("C",)].Pink_1
B_actin_control=grouped_trials[("C",)].B_actin
Abs_control=grouped_trials[("C",)].Absorbance

Pink_h2o2_100=grouped_trials[("100",)].Pink_1
B_actin_h2o2_100=grouped_trials[("100",)].B_actin
Abs_h2o2_100=grouped_trials[("100",)].Absorbance

Pink_h2o2_200=grouped_trials[("200",)].Pink_1
B_actin_h2o2_200=grouped_trials[("200",)].B_actin
Abs_h2o2_200=grouped_trials[("200",)].Absorbance

max_pink_values=[]
max_beta_actin_values=[]

push!(max_pink_values,maximum(Pink_control))
push!(max_pink_values,maximum(Pink_h2o2_100))
push!(max_pink_values,maximum(Pink_h2o2_200))

push!(max_beta_actin_values,maximum(B_actin_control))
push!(max_beta_actin_values,maximum(B_actin_h2o2_100))
push!(max_beta_actin_values,maximum(B_actin_h2o2_200))

triple_max_pink_values=max_pink_values.*3
triple_max_beta_actin_values=max_beta_actin_values.*3

std_pink_values=[]
std_beta_actin_values=[]

push!(std_pink_values,std(Pink_control))
push!(std_pink_values,std(Pink_h2o2_100))
push!(std_pink_values,std(Pink_h2o2_200))

push!(std_beta_actin_values,std(B_actin_control))
push!(std_beta_actin_values,std(B_actin_h2o2_100))
push!(std_beta_actin_values,std(B_actin_h2o2_200))

std_double_pink=std_pink_values.*2
std_double_beta_actin=std_beta_actin_values.*2

@model function normal_fit(data,index,triple_avg,double_std)
	μ ~ Di.Uniform(0,triple_avg[index])
	σ ~ Di.Uniform(0,double_std[index])
    data ~ Di.MvNormal(Fill(μ,length(data)),σ)
end

function distr_det(data,index,triple_avg,double_std)
	means=[]
	sigmas=[]
	model1 = normal_fit(data,index,triple_avg,double_std)
	chain = Turing.sample(model1,NUTS(0.65),1000)
	means=mean(chain[:μ])
	sigmas=mean(chain[:σ])
	#plot(chain)
	return means, sigmas
end

function get_mu_sigma(data,index,z_abs,tripple_avg,double_std)
	means, sigmas=distr_det(data,index,triple_avg,double_std)
	z=[-z_abs,z_abs]
	x=z.*sigmas[1].+means[1]
	plot(Di.Normal(means[1],sigmas[1]))
	pos_cutoff=	SF.erf(z[2]/sqrt(2)/2)
	neg_cutoff=1-pos_cutoff
	return x, means, sigmas
end

function create_list_of_outliers(data,index,z_abs,triple_avg,double_std)
	x, means, sigmas =get_mu_sigma(data,index,z_abs,triple_avg,double_std)
	outliers=[]
	for i=1:length(data)
		if (data[i]>x[2])||(data[i]<x[1])
			push!(outliers,data[i])
		end
	end
	return outliers
end

function MLE(data,index,triple_avg,double_std)
	means, sigmas = distr_det(data,index,triple_avg,double_std)
	likelyhoods=zeros(length(data))
	area_under_curve=zeros(length(data))
	for j=1:length(data)
		likelyhoods[j]= log(1/(sigmas[1]*sqrt(2*pi))*exp(-(data[j]-means[1])^2/(2*sigmas[1]^2)))
		area_under_curve[j]= (1/sqrt(2*pi*sigmas[1]^2)) * exp(-(means[1]-data[j]))
	end
	MLE_output=sum(likelyhoods)
	return MLE_output, likelyhoods, means
end

likleyhoods_after_below_threshold=[]
MLE_outputs_after_below_threshold=[]	

#=
	The function remove_data_below_threshold examines the likelihood of each data point in a set. If the likeleihood is below 5%, the data point will be removed and a new maximum likelihood estimation will be calculated. If there are multiple data points under this threshold then all of those ppints will be removed.

	The function remove_each_data_point_ind does the same thing as the previous function, but without any regard for a threshold. Instead this will help us see if there are any points within our dataset that don't qualify as "outliers" but still have a significiant impact on our data.  Note, this is evaluated one at a time so multiple data points will not be removed at once 
=#

function remove_data_below_threshold(data,index,triple_avg,double_std,threshold)
	MLE_output, likelyhoods, means= MLE(data,index,triple_avg,double_std)
	new_data=[]
	for i=1:length(likelyhoods)
		temp_data=[];
		temp_data=copy(data)
		MLE_output, likelyhoods, means= MLE(temp_data,index,triple_avg,double_std)
		exp_likelyhoods=exp.(likelyhoods)
		if exp_likelyhoods[i]<threshold#should be .05
			deleteat!(temp_data,i)
		end
		new_data=copy(temp_data)
	end
	triple_avg_new=zeros(index)
	double_std_new=zeros(index)
	triple_avg_new[index]=3*maximum(new_data)
	double_std_new[index]=2*std(new_data)
	MLE_output, likelyhoods, means= MLE(new_data,index,triple_avg_new,double_std_new)
	push!(likleyhoods_after_below_threshold, likelyhoods)
	push!(MLE_outputs_after_below_threshold, MLE_output)
end

#=
These results show us that removing the first data point gives us the highest maximum likelihood estimation. The higher the MLE, the better the fit of the data.Thus we qualify data point one from the 200uM experimental group as an outlier.
=#

function remove_each_data_point_ind(data,index,triple_avg,double_std)
	likleyhoods_after_each_data_is_removed=[]
	MLE_outputs_after_each_data_is_removed=[]
	MLE_output, likelyhoods, means= MLE(data,index,triple_avg,double_std)
	for i=1:length(likelyhoods)
		temp_data=[];
		temp_data=copy(data)
		MLE_output, likelyhoods, means= MLE(temp_data,index,triple_avg,double_std)
		deleteat!(temp_data,i)
		triple_avg_new=zeros(index)
		double_std_new=zeros(index)
		triple_avg_new[index]=3*maximum(temp_data)
		double_std_new[index]=2*std(temp_data)	
		MLE_output, likelyhoods, means= MLE(temp_data,index,triple_avg_new,double_std_new)
		push!(likleyhoods_after_each_data_is_removed, likelyhoods)
		push!(MLE_outputs_after_each_data_is_removed, MLE_output)
	end
	return likleyhoods_after_each_data_is_removed, MLE_outputs_after_each_data_is_removed
end

#MLE(B_actin_h2o2_200,3,triple_max_beta_actin_values,std_double_beta_actin)

#remove_data_below_threshold(B_actin_h2o2_200,3,triple_max_beta_actin_values,std_double_beta_actin,0.05)

#exp.(likleyhoods_after_below_threshold[1])

remove_each_data_point_ind(B_actin_h2o2_200,3,triple_max_beta_actin_values,std_double_beta_actin)

function get_final_data_set(b_actin_data,pink_data,index,triple_avg,double_std)
	likleyhoods_after_each_data_is_removed, MLE_outputs_after_each_data_is_removed =remove_each_data_point_ind(b_actin_data,index,triple_avg,double_std)
	mean_MLE=mean(MLE_outputs_after_each_data_is_removed);
	index_to_remove=[]
	new_b_actin_data=copy(b_actin_data)
	new_pink_data=copy(pink_data)
	for k=1:length(MLE_outputs_after_each_data_is_removed)
		if MLE_outputs_after_each_data_is_removed[k]>2*mean_MLE
			deleteat!(new_b_actin_data,k)
			deleteat!(new_pink_data,k)
		end
	end
	return new_b_actin_data,new_pink_data
end

#OneSampleTTest(Pink_control, new_pink_data)
#UnequalVarianceTTest(Pink_control, new_pink_data)

@model function normal_fit_no_index(data,triple_avg,double_std)
	μ ~ Di.Uniform(0,triple_avg)
	σ ~ Di.Uniform(0,double_std)
    data ~ Di.MvNormal(Fill(μ,length(data)),σ)
end

function distr_det_no_index(data,index,triple_avg,double_std)
	means=[]
	sigmas=[]
	#=if index==1
		avg = triple_avg[1]
		std = double_std[1]
	else =#
		max1 = 3* maximum(data)
		sd1 = 2*std(data)
	#end
	model1 = normal_fit_no_index(data,max1,sd1)
	chain = Turing.sample(model1,NUTS(0.65),1000)
	means=mean(chain[:μ])
	sigmas=mean(chain[:σ])
	#plot(chain)
	return means, sigmas
end

beta_c,pink_c= get_final_data_set(B_actin_control,Pink_control,1,triple_max_beta_actin_values,std_double_beta_actin)
means_c,sigma_c=distr_det_no_index(pink_c,1,triple_max_pink_values, std_double_pink)

beta_100,pink_100= get_final_data_set(B_actin_h2o2_100,Pink_h2o2_100,2,triple_max_beta_actin_values,std_double_beta_actin)
means_100,sigma_100=distr_det_no_index(pink_100,2,triple_max_pink_values, std_double_pink)

beta_200,pink_200= get_final_data_set(B_actin_h2o2_200,Pink_h2o2_200,3,triple_max_beta_actin_values,std_double_beta_actin)
means_200,sigma_200=distr_det_no_index(pink_200,3,triple_max_pink_values, std_double_pink)



