var virusTooltip = d3.tip()
	.direction('e')
	.attr('class', 'd3-tip')
	.offset([0, 12])
	.html(function(d) {
	
		string = "";
				
		// safe to assume the following attributes
		if (typeof d.strain != "undefined") {
			string += d.strain;
		}
		string += "<div class=\"smallspacer\"></div>";
		
		string += "<div class=\"smallnote\">";		
		
		// check if vaccine strain
		if (vaccineStrains.indexOf(d.strain) != -1) {
			string += "Vaccine strain<br>";
			var vaccine_date = new Date(vaccineChoice[d.strain]);

			string += "First chosen " + vaccine_date.toLocaleString("en-us", { month: "short" }) + " " + vaccine_date.getFullYear() + "<br>";
			string += "<div class=\"smallspacer\"></div>";
		}
		
		if (typeof d.region != "undefined") {
			string += d.region.replace(/([A-Z])/g, '$1');
		}
		if (typeof d.country != "undefined") {
			string += ", "+d.country.replace(/([A-Z])/g, '$1');
		}
		if (typeof d.date != "undefined") {
			string += ", " + d.date;
		}
		if (typeof d.host != "undefined") {
			string += "<br>Host: " + d.host;
		}
		if ((typeof d.accession != "undefined")) {
			string += "<br>Accession: " + d.accession;
		}
		if (typeof d.lab != "undefined") {
			if (d.lab != "") {
				string += "<br>Source: " + d.lab.substring(0,25);
				if (d.lab.length>25) string += '...';
			}
		}			
		string += "</div>";
		
		string += "<div class=\"smallspacer\"></div>";
				
		// following may or may not be present
		string += "<div class=\"smallnote\">";
		string+="Click to open genbank record" + "<br>";
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
		string+="Click to zoom into clade"
		if (typeof d.frequency != "undefined") {
			string += "<br>Frequency: " + (100 * d.frequency).toFixed(1) + "%"
		}
		if ((typeof d.aa_muts !="undefined")&&(d.aa_muts.length)){
			string+="<br>Mutations: "+d.aa_muts.replace(/,/g, ', ');
		}else if ((typeof d.nuc_muts !="undefined")&&(d.nuc_muts.length)){
			var tmp_muts = d.nuc_muts.split(',');
			var nmuts = tmp_muts.length;
			tmp_muts = tmp_muts.slice(0,Math.min(10, nmuts))
			string+="<br>Mutations: "+tmp_muts.join(', ');
			if (nmuts>10) {string+=' + '+ (nmuts-10) + ' more';}
		}
		return string;
	});
treeplot.call(linkTooltip);
