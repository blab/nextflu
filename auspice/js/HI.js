var predictedHI = true;
var correctVirus = true;
var correctPotency = true;
var focusNode;
/**
 * for each node, accumulate HI difference along branches
**/
function calcHIsubclade(node){
	node.HI_dist = node.parent.HI_dist+node.dHI;
	if (typeof node.children != "undefined") {
		for (var i=0; i<node.children.length; i++) {
		calcHIsubclade(node.children[i]);
		}
	}else{
		if (typeof node.avidity != "undefined" && correctVirus==false){
			node.HI_dist+=node.avidity;
		}
	}
};

function calcHIpred(node, rootNode){
	if (correctPotency){
		node.HI_dist = 0;
	}else{
		node.HI_dist=node.mean_potency;
	}
	if (typeof node.children != "undefined") {
		for (var i=0; i<node.children.length; i++) {
		calcHIsubclade(node.children[i]);
		}
	}
	var tmp_node = node;
	var pnode = tmp_node.parent;
	while (tmp_node.clade != rootNode.clade){
		pnode.HI_dist=tmp_node.HI_dist + tmp_node.dHI;
		if (typeof pnode.children != "undefined") {
			for (var i=0; i<pnode.children.length; i++) {
				if (tmp_node.clade!=pnode.children[i].clade){
					calcHIsubclade(pnode.children[i]);
				}
			}
		}
		tmp_node = pnode;
		pnode = tmp_node.parent;
	}
	if (correctVirus==false){
		node.HI_dist += node.avidity;
	}
};

function calcHImeasured(node, rootNode){
	console.log(node.strain+ ', mean_potency:'+node.mean_potency);
	for (var i=0; i<tips.length; i+=1){
		d = tips[i];
		if (typeof(node.HI_titers[d.clade])!="undefined"){
			d.HI_dist = node.HI_titers[d.clade]
			if (correctVirus){
				d.HI_dist -= d.avidity;
			}
			if (correctPotency){
				d.HI_dist -= node.mean_potency;
			}
		}else{
			d.HI_dist = 'NaN';
		}
	}
};

function tipHIvalid(d) {
	var vis = "visible";
	if (d.HI_dist =='NaN') {
		vis = "hidden";
	}
	return vis;
}

function getSera(tree_tips){
	return tree_tips.filter(function (d){return d.serum;})
}

d3.select("#serum")
	.on("change", colorByHIDistance);

d3.select("#virus")
	.on("change", colorByHIDistance);

d3.select("#HIPrediction")
	.on("change", colorByHIDistance);
