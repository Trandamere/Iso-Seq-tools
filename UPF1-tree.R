library(ape)
library(ggtree)
library(ggplot2)
library(RColorBrewer)
library(circlize)
library(tidytree)
library(ggtreeExtra)
library(ggstar)
library(treeio)
library(ggnewscale)
tree <- read.tree("upf1-tree")
p <- ggtree(tree, layout = "circular", branch.length = "none") + geom_tiplab(offset = 4.5, size = 2)+
#geom_text(aes(Label=node))
geom_highlight(node=96,fill="#3CB371")+geom_highlight(node=132,fill="#FF69B4")+geom_highlight(node=149,fill="#FFA500")+
geom_highlight(node=124,fill="#DDA0DD")
p

heatmap_df <- read.table("UPF1-olth-2.txt", header = TRUE, sep = "\t")
p3 <- p + geom_fruit(
  data = heatmap_df,
  geom=geom_tile,
  mapping=aes(x = Type2, y = species, fill = value),
  width = 1,
  height = 1,
  pwidth = 3,
  pangle = NA) + scale_fill_gradientn(
  colors = c("#00CED1", "#FFFACD", "#FF4500"),
  na.value = "white",
  name = "Similarity")
p3
p4<- p3+new_scale_fill()
df2 <- read.table("UPF1-Pfam-type-count.txt", header = TRUE, sep="\t")

p5=p4 + geom_fruit(
  data = df2,
  geom=geom_tile,
  mapping=aes(x = Type2, y = species, fill = category),
  width = 1,
  height = 1,
  pwidth = 5s) +
  scale_fill_manual(breaks=c("DUF6539-PF13086.8-PF13087.8", 
					"PF04851.17-PF09416.12-PF13086.8-PF13087.8-PF18141.3",
                              "PF05183.14-PF13086.8-PF13087.8",
					"PF09416.12-PF13086.8-PF13087.8-PF13604.8-PF18141.3",
					"PF09416.12-PF13086.8-PF13087.8-PF18141.3",
					"PF09416.12-PF13086.8-PF18141.3",
					"PF09416.12-PF18141.3",
					"PF12726.9-PF13086.8-PF13087.8",
					"PF12770.9-PF13374.8",
					"PF13086.8-PF13087.8",
					"PF13086.8-PF13087.8-PF18141.3",
					"PF13087.8",
					"PF13087.8-PF13604.8"), 
        values=c("#FF0000", "#00FF00", "#0000FF", "#FF00FF", 
		     "#FFFF00", "#00FFFF", "#FFA500", "#800080",
		     "#008000", "#000080", "#FF4500", "#00CED1", "#800000"))
p5












p5=gheatmap(p4, df2 , offset=1.5, width=0.1, font.size=3, 
        colnames_angle=45, hjust=0) +
    scale_fill_manual(breaks=c("DUF6539-PF13086.8-PF13087.8", 
					"PF04851.17-PF09416.12-PF13086.8-PF13087.8-PF18141.3",
                              "PF05183.14-PF13086.8-PF13087.8",
					"PF09416.12-PF13086.8-PF13087.8-PF13604.8-PF18141.3",
					"PF09416.12-PF13086.8-PF13087.8-PF18141.3",
					"PF09416.12-PF13086.8-PF18141.3",
					"PF09416.12-PF18141.3",
					"PF12726.9-PF13086.8-PF13087.8",
					"PF12770.9-PF13374.8",
					"PF13086.8-PF13087.8",
					"PF13086.8-PF13087.8-PF18141.3",
					"PF13087.8",
					"PF13087.8-PF13604.8"), 
        values=c("#FF0000", "#00FF00", "#0000FF", "#FF00FF", 
		     "#FFFF00", "#00FFFF", "#FFA500", "#800080",
		     "#008000", "#000080", "#FF4500", "#00CED1", "#800000"))
p5












p5 <- gheatmap(p4, df2, offset=15, width=.3,
         colnames_angle=90, colnames_offset_y = .25) +
    scale_fill_viridis_d(option="D", name="discrete\nvalue")
p5



#第三圈
df2 <- read.table("UPF1-Pfam-type-count.txt", header = TRUE, sep="\t")
rownames(df2) <- tree$tip.label
colnames(df2) <- df2$Sites

p2=gheatmap(p, df2 , offset=0.5, width=0.1, font.size=3, 
        colnames_angle=-45, hjust=0) +
    scale_fill_manual(breaks=c("DUF6539-PF13086.8-PF13087.8", 
					"PF04851.17-PF09416.12-PF13086.8-PF13087.8-PF18141.3",
                              "PF05183.14-PF13086.8-PF13087.8",
					"PF09416.12-PF13086.8-PF13087.8-PF13604.8-PF18141.3",
					"PF09416.12-PF13086.8-PF13087.8-PF18141.3",
					"PF09416.12-PF13086.8-PF18141.3",
					"PF09416.12-PF18141.3",
					"PF12726.9-PF13086.8-PF13087.8",
					"PF12770.9-PF13374.8",
					"PF13086.8-PF13087.8",
					"PF13086.8-PF13087.8-PF18141.3",
					"PF13087.8",
					"PF13087.8-PF13604.8"), 
        values=c("#FF0000", "#00FF00", "#0000FF", "#FF00FF", 
		     "#FFFF00", "#00FFFF", "#FFA500", "#800080",
		     "#008000", "#000080", "#FF4500", "#00CED1", "#800000"))
p2
# 在圆形发育树上添加热图




scale_fill_manual(
  values = c("DUF6539-PF13086.8-PF13087.8"="#FF0000","PF04851.17-PF09416.12-PF13086.8-PF13087.8-PF18141.3"="#00FF00",
  "PF05183.14-PF13086.8-PF13087.8"="#0000FF","PF09416.12-PF13086.8-PF13087.8-PF13604.8-PF18141.3"="#FF00FF",
  "PF09416.12-PF13086.8-PF13087.8-PF18141.3"="#FFFF00",
  "PF09416.12-PF13086.8-PF18141.3"="#00FFFF","PF09416.12-PF18141.3"="#FFA500","PF12726.9-PF13086.8-PF13087.8"="#800080",
  "PF12770.9-PF13374.8"="#008000","PF13086.8-PF13087.8"="#000080","PF13086.8-PF13087.8-PF18141.3"="#FF4500",
  "PF13087.8"="#00CED1","PF13087.8-PF13604.8"="#800000")
