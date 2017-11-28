using HDF5, JLD, DataFrames, Distributions, MLBase, StatsBase
using DecisionTree: mean_squared_error

function whos_user(m::Module=current_module())
    for v in sort(names(m))
        s = string(v)
        if isdefined(m, v) && summary(eval(m, v)) != "Module" && s != "whos_user"
            println(s)
        end
    end
end

function text_to_numb(A)
    novo = Array{Float64}(size(A)[1],1);
    n = unique(dropna(A));
    for i in 1:size(A)[1]
        if isna(A[i])
            novo[i] = -1
        else
            power = findfirst(n,A[i]) - 1; 
            novo[i] = 2^power - 1;
        end
    end
    return novo;
end

function break_this_baby(A)
	lefty = Array{Float64}(size(A)[1],1);
	righty = Array{Float64}(size(A)[1],1);
	for i in 1:size(A)[1]
		n = search(A[i], "_")[1] - 1;
		lefty[i] = parse(Float64, A[i][2:n]);
		m = search(A[i], "]")[1] - 1;
		n = n+2;
		righty[i] = parse(Float64, A[i][n:m]);
	end
	return lefty, righty;
end

function break_this_weirdo(A,k)
    out = zeros(size(A)[1],k);
    for i in 1:size(A)[1]
        if !isna(A[i])
            reg = matchall(r"[0-9]+\.[0-9]+", A[i])
            for j in 1:min(k,length(reg))
                out[i,j] = parse(Float64, reg[j])
            end
        end
    end
    return out;
end

function root_mean_squared_error(A,B)
  return sqrt(mean_squared_error(A,B))
end

function confusion_matrix_regression(actual,predicted,k)
	@assert length(actual) == length(predicted)
	out = round(Int64,zeros(k,k))
	frame = DataFrame();
	names = Array{ASCIIString}(1,k)
	ranger = collect(linspace(min(minimum(actual),minimum(predicted)),max(maximum(actual),maximum(predicted)),k+1))
	for i in 1:length(actual)
		j = min(searchsortedlast(ranger, actual[i]),k)
		l = min(searchsortedlast(ranger, predicted[i]),k)
		out[j,l] = out[j,l] + 1
	end
	
	for i in 1:k
       names[i] = string(trunc(ranger[i],2))*"-"*string(trunc(ranger[i+1],2))
       frame[symbol(names[i])] = out[:,i]
    end
    frame[:TOTAL] = convert(DataArray,sum(out,2))[:]
	return frame
end

function cross_validate_forests(data, Y, k, num_of_trees, num_of_features)
	p = randperm(size(data)[1])
	data = data[p,:]
	Y = Y[p,:]
	ranger = round(Int64,floor(collect(linspace(1,size(data)[1],k+1))))
	RMSEs = Array{Float64}(k,1)
	BestModel = 0
	for i in 1:k
		train = trues(size(data)[1])
		train[ranger[i]:ranger[i+1]] = false

		println("\nRUN #",i)
		display(Dates.format(now(), "e, dd u yyyy HH:MM"))
		@time Model = build_forest(Y[train],data[train,:],num_of_features,num_of_trees);
		Prediction = apply_forest(Model, data[!train,:]);
		RMSEs[i] = root_mean_squared_error(Y[!train], Prediction);
		if(i == 1)
			BestModel = (deepcopy(Model), RMSEs[i])
		elseif RMSEs[i] < BestModel[2]
			BestModel = (deepcopy(Model), RMSEs[i])
		end
		println("\n", Model)
		println("Features: ",num_of_features)
		println("RMSE: ",RMSEs[i])
		display(confusion_matrix_regression(Y[!train],Prediction,10))
	end
	println("\nMean RMSE: ",mean(RMSEs))
	println("Min RMSE: ",minimum(RMSEs))
	println("Max RMSE: ",maximum(RMSEs))

	return RMSEs, BestModel
end

function splitbit_tr_te(bita, split)
	tr = falses(length(bita))
	te = falses(length(bita)) 
	for i in 1:length(bita)
		if bita[i] == true
			if rand(Bernoulli(split)) == 1
				tr[i] = true
			else
				te[i] = true
			end
		end
	end
	return tr, te
end

function splitbit_tr_te_xval(bita, k)
	tr = falses(length(bita),k)
	te = falses(length(bita),k)
	split = floor(sum(bita)/k)
	break_here = 0;
	for j in 1:k
		count = 0;
		for i in 1:length(bita)
			if bita[i] == true
				count = count + 1;
				if count > break_here && count <= j*split
					te[i,j] = true
				else
					tr[i,j] = true
				end
			end
		end
		break_here = break_here + split;
	end
	return tr, te
end

const NO_BEST=(0,0)
function rank_parameters{T<:Float64, U<:Real}(labels::Vector{T}, features::Matrix{U}, numoffeatures::Int)
    nr, nf = size(features)
    the_list = []
    ignore = []
    value = zeros(nf)
    for l in 1:numoffeatures
        best = NO_BEST
        best_val = -Inf
        for i in 1:nf
            if !(i in ignore)
                if nr > 100
                    features_i = features[:,i]
                    domain_i = quantile(features_i, linspace(0.01, 0.99, 99))
                    labels_i = labels
                else
                    ord = sortperm(features[:,i])
                    features_i = features[ord,i]
                    domain_i = features_i
                    labels_i = labels[ord]
                end

                for thresh in domain_i[2:end]
                    value = DecisionTree._mse_loss(labels_i, features_i, thresh)
                    if value > best_val
                        best_val = value
                        best = (i, value, thresh)
                    end
                end
            end
        end
        push!(the_list,best)
        push!(ignore,best[1])
    end
    return the_list
end

function translate_to_names(bestF)
    output = AbstractString[]
    for a in bestF
        push!(output, fullDataNames[a[1]])
    end
    return output
end

function unique_indices(bestP)
    index = falses(length(bestP))
    for i in unique(bestP)
        for j in 1:length(bestP)
            if bestP[j] == i
                index[j] = !index[j]
                break
            end
        end
    end
    return index
end

function accuracy(predictions, classes)
    @assert length(predictions) == length(classes)
    class = sort(unique(classes))
    right = 0
    for i = 1:size(predictions)[1]
        predictedClass = menor = Inf
        for j in class
            if abs(j - predictions[i]) < menor
                menor = abs(j - predictions[i])
                predictedClass = j
            end
        end
        if predictedClass == classes[i]
            right = right+1
        end
    end
    return (right/length(predictions))*100
end

function get_weight(spirality)
    #weight = 2*(spirality.^2)-2.*spirality+1 #parabola 
    weight = (0.5.*spirality).+0.5
    return weight
end