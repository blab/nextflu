var HImodel = 'measured';
var correctVirus = true;
var correctPotency = true;
var focusNode;
var activeSera = {};
/**
 * for each node, accumulate HI difference along branches
**/
function calcHIsubclade(node){
	node.HI_dist_tree = node.parent.HI_dist_tree+node.attr.dTiter;
	if (typeof node.children != "undefined") {
		for (var i=0; i<node.children.length; i++) {
		calcHIsubclade(node.children[i]);
		}
	}else{
		if (typeof titer_tree_model["avidity"][node.strain] != "undefined" && correctVirus==false){
			node.HI_dist_tree+=titer_tree_model["avidity"][node.strain];
		}
	}
};

function calcHItree(node, rootNode){
	if (correctPotency){
		node.HI_dist_tree = 0;
	}else{
		node.HI_dist_tree=titer_tree_model["potency"][node.strain].mean_potency;
	}
	if (typeof node.children != "undefined") {
		for (var i=0; i<node.children.length; i++) {
		calcHIsubclade(node.children[i]);
		}
	}
	var tmp_node = node;
	var pnode = tmp_node.parent;
	while (tmp_node.strain != rootNode.strain){
		pnode.HI_dist_tree=tmp_node.HI_dist_tree + tmp_node.attr.dTiter;
		if (typeof pnode.children != "undefined") {
			for (var i=0; i<pnode.children.length; i++) {
				if (tmp_node.strain!=pnode.children[i].strain){
					calcHIsubclade(pnode.children[i]);
				}
			}
		}
		tmp_node = pnode;
		pnode = tmp_node.parent;
	}
	if (correctVirus==false){
		node.HI_dist_tree += titer_tree_model["avidity"][node.strain];
	}
};

function calcHImeasured(node, rootNode){
	console.log(node.strain+ ', mean_potency: '+titer_subs_model["potency"][node.strain].mean_potency);
	console.log("correcting for virus effect: "+correctVirus);
	console.log("correction for serum effect: "+correctPotency);
	var tmptt;
	for (var i=0; i<tips.length; i+=1){
		d = tips[i];
		if (typeof HI_titers[node.strain][d.strain] != "undefined"){
			var tmptt = HI_titers[node.strain][d.strain];
			var tmp_HI=0;
			var serum_count=0;
			for (var tmp_serum in tmptt){
				if (activeSera[tmp_serum]){
					if (correctPotency&&(d.strain!=focusNode.strain)){
						tmp_HI += tmptt[tmp_serum][0]-titer_subs_model["potency"][node.strain].mean_potency
					}else{
						tmp_HI += tmptt[tmp_serum][0];
					}
					serum_count+=1;
				}
			}
			if (serum_count){
				d.HI_dist_meas = tmp_HI/serum_count
				if (correctVirus){
					d.HI_dist_meas -= titer_subs_model["avidity"][d.strain];
				}
			}else{
				d.HI_dist_meas = 'NaN';
			}
		}else{
			d.HI_dist_meas = 'NaN';
		}
	}
};

function get_mutations(node1, node2){
	var gt1, gt2,muts=[];
	for (var gene in cladeToSeq[node1.strain]){
		var gene_length = cladeToSeq["root"][gene].length;
		if (gene!='nuc'){
			for (var pos=0; pos<gene_length; pos+=1){
				gt1 = stateAtPosition(node1.strain, gene, pos)
				gt2 = stateAtPosition(node2.strain, gene, pos)
				if (gt1!=gt2){
					muts.push(gene+':'+gt1+(pos+1)+gt2);
				}
			}
		}
	}
	return muts;
}

function calcHImutations(node){
	console.log(node.strain+ ', mean_potency:'+titer_subs_model["potency"][node.strain].mean_potency);
	nodes.map(function(d){
		var mutations = get_mutations(node, d);
		if (correctPotency){
			d.HI_dist_mut=0;
		}else{
			d.HI_dist_mut=titer_subs_model["potency"][node.strain].mean_potency;
		}
		for (var mi=0; mi<=mutations.length; mi++){
			var mut = mutations[mi];
			if ((typeof mut != "undefined")&&(typeof titer_subs_model["substitution"][mut]!="undefined")){
				d.HI_dist_mut += titer_subs_model["substitution"][mut];
			}
		}
		if ((correctVirus==false)&&(typeof d.avidity != "undefined")){
			d.HI_dist_mut += titer_subs_model["avidity"][d.strain];
		}
	});
};

function getSera(tree_tips){
	return tree_tips.filter(function (d){return typeof(HI_titers[d.strain]) != "undefined";})
}

d3.select("#serum")
	.on("change", colorByHIDistance);

d3.select("#virus")
	.on("change", colorByHIDistance);

d3.select("#HImodel_measured")
	.on("click", colorByHIDistance);
d3.select("#HImodel_mutation")
	.on("click", colorByHIDistance);
d3.select("#HImodel_tree")
	.on("click", colorByHIDistance);

var HI_titers, titer_tree_model, titer_subs_model;
var structure_HI_mutations;
if (useTiters) {
	d3.json(path + file_prefix + "titers.json", function(error, json){
		HI_titers = json;
	});
	d3.json(path + file_prefix + "titer-tree-model.json", function(error, json){
		titer_tree_model = json;
	});
	d3.json(path + file_prefix + "titer-sub-model.json", function(error, json){
		titer_subs_model = json;
	});
}
