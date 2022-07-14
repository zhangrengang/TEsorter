#!/bin/env Rscript
args = commandArgs(T)
treefile = args[1]
mapfile = args[2]
outfig = args[3]
if (is.na(outfig)) {outfig = paste0(treefile, '.png')}

branch_color = 'Clade'
library(ggplot2)
library(ggtree)
library(treeio)

map = read.table(mapfile, head=T, fill=T, comment.char='!', sep="\t")
tree <- read.tree(file = treefile)

split_id <- function(x) {
	x = strsplit(x, '#')[[1]][1]
	return(x)
}
format_id <- function(x1, x2, x3, x4) {
	x1 = sapply(x1, split_id)	
	x1 = gsub('\\W+', '_', x1)
	x = paste(x1, x2, x3, x4, sep='_')
	return(x)
}

if (branch_color == 'Clade') {
	clades = sort(unique(map$Clade))
	head(tree)
	tree$tip.label = gsub('\\W+', '_', tree$tip.label)
	print(clades)
	grp = list()
	for (clade in clades){
			if (clade=='unknown') { next }
			labels = map[which(map$Clade==clade), ]
			clade = paste(labels$Superfamily, labels$Clade, sep='/')[1]
			labels = format_id(labels$X.TE, labels$Order, labels$Superfamily, labels$Clade)
			if (! any(labels %in% tree$tip.label)) {next}
	#		print(clade)
			grp[[clade]] = labels
	}
	clades = sort(names(grp))
	tree3 = groupOTU(tree, grp, 'Clade')
	p = ggtree(tree3 , aes(color=Clade) ,  layout='circular' ) +
	  theme(legend.position="right")  + 
	  scale_fill_manual(values=c('#f9c00c','#00b9f1','#7200da','#f9320c','#980000','#00ffff','#0000ff','#ff0000','#4a86e8','#ff9900','#ffff00','#00ff00','#9900ff','#ff00ff','#20124d','#274e13','#000000','#cccccc','#7f6000','#a64d79','#6aa84f','#fff2cc','#47a952','#3ea6b6','#a5b805','#8f9276','#ca8d7c')) + scale_colour_discrete(limits=clades, labels=clades) +
	  guides(colour=guide_legend(order = 1), fill=guide_legend(order = 2))

} else {	# branch_color == 'Taxon'
	taxa = sort(unique(map$Taxon))

	grp = list()
	for (taxon in taxa){
			labels = map[which(map$Taxon==taxon), ]
			labels = labels$label
			grp[[taxon]] = labels
	}
	tree3 = groupOTU(tree, grp, 'Taxon')
	map3 = data.frame(label=map$label, Clade=map$Clade)
	p = ggtree(tree3 , aes(color=Taxon) , branch.length='none', layout='circular' ) %<+% map3 +
	  theme(legend.position="right")  + 
	  scale_colour_manual(values=c('#f9c00c','#00b9f1','#7200da','#f9320c','#980000','#00ffff','#0000ff','#ff0000','#4a86e8','#ff9900','#ffff00','#00ff00','#9900ff','#ff00ff','#20124d','#274e13','#000000','#cccccc','#7f6000','#a64d79','#6aa84f','#fff2cc','#47a952','#3ea6b6','#a5b805','#8f9276','#ca8d7c'),limits=taxa, labels=taxa) +
	  geom_tippoint(aes(fill=Clade), pch=21, stroke=0, size=1.2, color='#00000000') +
	  scale_fill_hue(l=35) +
	  guides(colour=guide_legend(order = 1), fill=guide_legend(order = 2))

}
position = c(1.35,0.9)
p = p + theme(plot.margin=margin(0,0,0,0)) +
	theme(legend.position=position, legend.justification=position) +
	theme(legend.background=element_blank(), legend.key=element_blank()) +
	theme(legend.text=element_text(size=18), legend.title=element_text(size=20))


ggsave(outfig, p, width=13.5, height=8.4, dpi=300, units="in")
