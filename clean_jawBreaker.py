from pymol import cmd, stored

def cleanSim(d_resid, d_color, e_resid, e_color):
    """
    Cleans the simulation, colors drug d_color, 
    excipient e_color, and sets Zoom
    
    @param d_resid: input 3 char residue ID of drug
    @param d_color: input 6 char color code hash for drug
    @param e_resid: input 3 char residue ID of excipient
    @param e_color: 6 char color code hash for excipient

    
    Returns
    -------
    None.

    """
 
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", "on")
    cmd.set("orthoscopic", "1")
    cmd.set("depth_cue", "0")
    cmd.remove("solvent")
    cmd.show(representation = "spheres", selection = "all")
    cmd.color(f"{d_color}",  f"resname {d_resid}")
    cmd.color(f"{e_color}",  f"resname {e_resid}")
    cmd.color("yellow", "resname STK")

    cmd.center(f"resname {d_resid}")
    cmd.orient()
    cmd.zoom("center", "65")
 

cmd.extend("cleanSim", cleanSim)