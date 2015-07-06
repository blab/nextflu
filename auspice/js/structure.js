var myapplett;

$(document).ready(
	function() {
	console.log('drawing structure');
	Info = {
		width: 500,
		height: 500,
		debug: false,
		j2sPath: "/js/j2s",
		color: "white",
		disableJ2SLoadMonitor: true,
		disableInitialConsole: true,
		addSelectionOptions: false,
		use: "HTML5",
		readyFunction: null,
		script: "load /data/5HMG.pdb; cpk off; wireframe off; trace;zoom on; zoom 115;set showhydrogens off; color background white; select ligand; trace off; spin on; set frank off; set echo bottom left; color echo gray; font echo 14 arial;select 154:a,154:c,154:b;spacefill 200;color orange; select 239:a,239:c,239:e;spacefill 200;color red;"
	}

	myapplett = $("#HA_struct").html(Jmol.getAppletHtml("jmolApplet0",Info));

	$("#HA_buttons").html(
		Jmol.jmolButton(jmolApplet0, "spin on","Spin ON")
	 +Jmol.jmolButton(jmolApplet0, "spin off","Spin OFF")
	 +Jmol.jmolButton(jmolApplet0, "write PNGJ flusurver.png","Save IMAGE")
	 );
	console.log("script select 159;  spacefill 200;  color blue;");
	Jmol.script(myapplett, "script select 159;  spacefill 200;  color blue;");
});
