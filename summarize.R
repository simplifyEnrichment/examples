
library(GetoptLong)
library(simplifyGO)

if(normalizePath("~") == "/Users/jokergoo") {
	setwd("/Users/jokergoo/project/simplifyGO_figures")
} else {
    setwd("/icgc/dkfzlsdf/analysis/B080/guz/repo/simplifyGO_figures")
}

concordance_mat_lt = list()
diff_score_df = data.frame(method = character(0), value = numeric(0))
cluster_number_df = data.frame(method = character(0), value = numeric(0))
cluster_number2_df = data.frame(method = character(0), value = numeric(0))
block_mean_df = data.frame(method = character(0), value = numeric(0))
for(i in 1:100) {
	lt = readRDS(qq("random_BP/rds/clt_random_BP_@{i}.rds"))

	concordance_mat_lt[[i]] = simplifyGO:::compare_methods_calc_concordance(lt[[2]])
	
	x = sapply(lt[[2]], function(x) difference_score(lt[[1]], x))
	diff_score_df = rbind(diff_score_df, data.frame(method = names(x), value = x))
	
	x = sapply(lt[[2]], function(x) length(table(x)))
	cluster_number_df = rbind(cluster_number_df, data.frame(method = names(x), value = x))
	
	x = sapply(lt[[2]], function(x) {tb = table(x); sum(tb >= 5)})
	cluster_number2_df = rbind(cluster_number2_df, data.frame(method = names(x), value = x))
	
	x = sapply(lt[[2]], function(x) simplifyGO:::block_mean(lt[[1]], x))
	block_mean_df = rbind(block_mean_df, data.frame(method = names(x), value = x))
}

diff_score_df[, 1] = factor(diff_score_df[, 1], levels = simplifyGO:::ALL_CLUSTERING_METHODS)
cluster_number_df[, 1] = factor(cluster_number_df[, 1], levels = simplifyGO:::ALL_CLUSTERING_METHODS)
cluster_number2_df[, 1] = factor(cluster_number2_df[, 1], levels = simplifyGO:::ALL_CLUSTERING_METHODS)
block_mean_df[, 1] = factor(block_mean_df[, 1], levels = simplifyGO:::ALL_CLUSTERING_METHODS)

library(ggplot2)
library(grid)
p1 = ggplot(diff_score_df, aes(x = method, y = value)) +
	geom_boxplot() + ylab("Difference score") +
	theme(axis.title.x = element_blank(), axis.text.x = element_blank())

p2 = ggplot(cluster_number_df, aes(x = method, y = value)) +
	geom_boxplot() + ylab("Mean cluster size") +
	theme(axis.title.x = element_blank(), axis.text.x = element_blank())

p3 = ggplot(block_mean_df, aes(x = method, y = value)) +
	geom_boxplot() + ylab("Block mean") +
	theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

library(svglite)
# png("random_BP_boxplot.png", width = 400, height = 600, res = 100)
svglite("random_BP_boxplot.svg", width = 4, height = 6)
grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), ggplotGrob(p3), size = "last"))
dev.off()


arr = array(dim = c(dim(concordance_mat_lt[[1]]), length(concordance_mat_lt)))
dimnames(arr) = c(dimnames(concordance_mat_lt[[1]]), list(1:length(concordance_mat_lt)))
for(i in 1:length(concordance_mat_lt)) {
	arr[, , i] = concordance_mat_lt[[i]]
}

m_max = apply(arr, 1:2, max)
m_min = apply(arr, 1:2, min)
m_mean = apply(arr, 1:2, mean)

library(ComplexHeatmap)
# png("random_BP_similarity.png", width = 500, height = 500)
svglite("random_BP_similarity.svg", width = 5, height = 5)
Heatmap(m_mean, name = "similarity", cell_fun = function(j, i, x, y, w, h, f) {
	grid.rect(x, y, w, h, gp = gpar(col = "white", fill = "#EFEFEF"))
	grid.circle(x, y, r = min(w, h)*0.9*0.5*m_max[i, j], gp = gpar(col = f, fill = f))
	grid.circle(x, y, r = min(w, h)*0.9*0.5*m_mean[i, j], gp = gpar(fill = NA))
	grid.circle(x, y, r = min(w, h)*0.9*0.5*m_min[i, j], gp = gpar(col = f, fill = "#EFEFEF"))
}, rect_gp = gpar(type = "none"))
dev.off()

saveRDS(list(diff_score_df = diff_score_df,
	         cluster_number_df = cluster_number_df,
	         cluster_number2_df = cluster_number2_df,
	         block_mean_df = block_mean_df,
	         arr = arr),
        file = "random_BP_results.rds",
        compress = "xz")

file.copy("random_BP_results.rds", "../simplifyGO/vignettes/random_BP_results.rds", overwrite = TRUE)

##################################
bp_list = readRDS("bp_list.rds")

go_list = lapply(bp_list, function(x) {
	names(x)[p.adjust(x, "BH") < 0.05]
})

n = sapply(go_list, length)
go_list = go_list[n > 100 & n < 1000]

nm = names(go_list)
j = 0
concordance_mat_lt = list()
diff_score_df = data.frame(method = character(0), value = numeric(0))
cluster_number_df = data.frame(method = character(0), value = numeric(0))
cluster_number2_df = data.frame(method = character(0), value = numeric(0))
block_mean_df = data.frame(method = character(0), value = numeric(0))

for(i in seq_along(go_list)) {
	oe = try(lt <- readRDS(qq("EBI_Expression_Atlas/rds/clt_@{nm[i]}.rds")))

	if(inherits(oe, "try-error")) {
		next
	}

	j = j + 1
	concordance_mat_lt[[j]] = simplifyGO:::compare_methods_calc_concordance(lt[[2]])
	
	x = sapply(lt[[2]], function(x) difference_score(lt[[1]], x))
	diff_score_df = rbind(diff_score_df, data.frame(method = names(x), value = x))
	
	x = sapply(lt[[2]], function(x) length(table(x)))
	cluster_number_df = rbind(cluster_number_df, data.frame(method = names(x), value = x))
	
	x = sapply(lt[[2]], function(x) {tb = table(x); sum(tb >= 5)})
	cluster_number2_df = rbind(cluster_number2_df, data.frame(method = names(x), value = x))
	
	x = sapply(lt[[2]], function(x) simplifyGO:::block_mean(lt[[1]], x))
	block_mean_df = rbind(block_mean_df, data.frame(method = names(x), value = x))
}

diff_score_df[, 1] = factor(diff_score_df[, 1], levels = simplifyGO:::ALL_CLUSTERING_METHODS)
cluster_number_df[, 1] = factor(cluster_number_df[, 1], levels = simplifyGO:::ALL_CLUSTERING_METHODS)
cluster_number2_df[, 1] = factor(cluster_number2_df[, 1], levels = simplifyGO:::ALL_CLUSTERING_METHODS)
block_mean_df[, 1] = factor(block_mean_df[, 1], levels = simplifyGO:::ALL_CLUSTERING_METHODS)

library(ggplot2)
p1 = ggplot(diff_score_df, aes(x = method, y = value)) +
	geom_boxplot() + ylab("Difference score") +
	theme(axis.title.x = element_blank(), axis.text.x = element_blank())

p2 = ggplot(cluster_number_df, aes(x = method, y = value)) +
	geom_boxplot() + ylab("Mean cluster size") +
	theme(axis.title.x = element_blank(), axis.text.x = element_blank())

p3 = ggplot(block_mean_df, aes(x = method, y = value)) +
	geom_boxplot() + ylab("Block mean") +
	theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

svglite("EBI_Expression_Atlas_boxplot.svg", width = 4, height = 6)
# png("EBI_Expression_Atlas_boxplot.png", width = 400, height = 600, res = 100)
grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), ggplotGrob(p3), size = "last"))
dev.off()

arr = array(dim = c(dim(concordance_mat_lt[[1]]), length(concordance_mat_lt)))
dimnames(arr) = c(dimnames(concordance_mat_lt[[1]]), list(1:length(concordance_mat_lt)))
for(i in 1:length(concordance_mat_lt)) {
	arr[, , i] = concordance_mat_lt[[i]]
}

m_max = apply(arr, 1:2, max)
m_min = apply(arr, 1:2, min)
m_mean = apply(arr, 1:2, mean)

svglite("EBI_Expression_Atlas_similarity.svg", width = 5, height = 5)
# png("EBI_Expression_Atlas_similarity.png", width = 500, height = 500)
Heatmap(m_mean, cell_fun = function(j, i, x, y, w, h, f) {
	grid.rect(x, y, w, h, gp = gpar(col = "white", fill = "#EFEFEF"))
	grid.circle(x, y, r = min(w, h)*0.9*0.5*m_max[i, j], gp = gpar(col = f, fill = f))
	grid.circle(x, y, r = min(w, h)*0.9*0.5*m_mean[i, j], gp = gpar(fill = NA))
	grid.circle(x, y, r = min(w, h)*0.9*0.5*m_min[i, j], gp = gpar(col = f, fill = "#EFEFEF"))
}, rect_gp = gpar(type = "none"))
dev.off()

saveRDS(list(diff_score_df = diff_score_df,
	         cluster_number_df = cluster_number_df,
	         cluster_number2_df = cluster_number2_df,
	         block_mean_df = block_mean_df,
	         arr = arr),
        file = "EBI_Expression_Atlas_results.rds",
        compress = "xz")

file.copy("EBI_Expression_Atlas_results.rds", "../simplifyGO/vignettes/EBI_Expression_Atlas_results.rds", overwrite = TRUE)
