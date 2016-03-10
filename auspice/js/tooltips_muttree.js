var virusTooltip = d3.tip()
	.direction('e')
	.attr('class', 'd3-tip')
	.offset([0, 12])
	.html(function(d) {

		string = "";
		if (typeof vaccineChoice[d.strain] != "undefined"){
			string+="REFERENCE STRAIN<br>";
		}

		// safe to assume the following attributes
		if (typeof d.strain != "undefined") {
			string += d.strain;
		}
		string += "<div class=\"smallspacer\"></div>";

		string += "<div class=\"smallnote\">";

		if (typeof d.country != "undefined") {
			string += d.country.replace(/([A-Z])/g, ' $1');
		}
		if (typeof d.date != "undefined") {
			string += ", " + d.date;
		}
		if (typeof d.isolate_id != "undefined") {
			string += "<br>Isolate: " + d.isolate_id;
		}
		if (typeof d.accession != "undefined") {
			string += "<br>GISAID ID: EPI" + d.accession;
		}
		if ((typeof d.subtype != "undefined")&&(d.subtype != "Unknown")) {
			if (d.orig_lab != "") {
				string += "<br>Subtype: " + d.subtype;
			}
		}
		if ((typeof d.host != "undefined")&&(d.host != "Unknown")) {
			if (d.orig_lab != "") {
				string += "<br>host group: " + d.host;
			}
		}
		if ((typeof d.passage != "undefined")&&(d.passage != "Unknown")) {
			if (d.orig_lab != "") {
				string += "<br>passage: " + d.passage;
			}
		}
		string += "</div>";

		string += "<div class=\"smallspacer\"></div>";

		// following may or may not be present
		string += "<div class=\"smallnote\">";
		if (typeof d.ep != "undefined") {
			string += "Epitope distance: " + d.ep + "<br>";
		}
		if (typeof d.ne != "undefined") {
			string += "Non-epitope distance: " + d.ne + "<br>";
		}
		if (typeof d.rb != "undefined") {
			string += "Receptor binding distance: " + d.rb + "<br>";
		}
		if (typeof d.LBI != "undefined") {
			string += "Local branching index: " + d.LBI.toFixed(3) + "<br>";
		}
		string += "</div>";
		return string;
	});
treeplot.call(virusTooltip);

var linkTooltip = d3.tip()
	.direction('e')
	.attr('class', 'd3-tip')
	.offset([0, 12])
	.html(function(d) {
		string = ""
		if (typeof d.frequency != "undefined") {
			string += "Frequency: " + (100 * d.frequency).toFixed(1) + "%"
		}
		if ((typeof d.aa_muts !="undefined")&&(d.aa_muts.length)){
			string+="<br>Mutations: "+d.aa_muts.replace(/,/g, ', ');
		}else if ((typeof d.nuc_muts !="undefined")&&(d.nuc_muts.length)){
			string+="<br>Mutations: "+d.nuc_muts.replace(/,/g, ', ');
		}
		return string;
	});
treeplot.call(linkTooltip);
