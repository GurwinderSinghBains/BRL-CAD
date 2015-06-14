/*                         G - X X X . C
 * BRL-CAD
 *
 * Copyright (c) 1993-2014 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this file; see the file named COPYING for more
 * information.
 *
 */
/** @file conv/g-xxx.c
 *
 * Sample code for converting BRL-CAD models to some other format.
 * This code assumes that your receiving format can handle CSG
 * primitives and Boolean trees with transformation matrices
 *
 */

#include "common.h"

/* system headers */
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "bio.h"

/* interface headers */
#include "vmath.h"
#include "bu/getopt.h"
#include "nmg.h"
#include "rt/geom.h"
#include "raytrace.h"
#include "wdb.h"


struct user_data {
    long int data;
    struct bn_tol tol;
};

/*

torus macro

*/


void torusmacro();


/* This routine just produces an ascii description of the Boolean tree.
 * In a real converter, this would output the tree in the desired format.
 */
void
describe_tree(union tree *tree,
	      struct bu_vls *str)
{
    struct bu_vls left = BU_VLS_INIT_ZERO;
    struct bu_vls right = BU_VLS_INIT_ZERO;
    const char op_xor='^';
    char op='\0';

    BU_CK_VLS(str);

    if (!tree) {
	/* this tree has no members */
	bu_vls_strcat(str, "-empty-");
	return;
    }

    RT_CK_TREE(tree);

    /* Handle all the possible node types.
     * the first four are the most common types, and are typically
     * the only ones found in a BRL-CAD database.
     */
    switch (tree->tr_op) {
	case OP_DB_LEAF:	/* leaf node, this is a member */
	    /* Note: tree->tr_l.tl_mat is a pointer to a
	     * transformation matrix to apply to this member
	     */
	    bu_vls_strcat(str,  tree->tr_l.tl_name);
	    break;
	case OP_UNION:		/* union operator node */
	    op = DB_OP_UNION;
	    goto binary;
	case OP_INTERSECT:	/* intersection operator node */
	    op = DB_OP_INTERSECT;
	    goto binary;
	case OP_SUBTRACT:	/* subtraction operator node */
	    op = DB_OP_SUBTRACT;
	    goto binary;
	case OP_XOR:		/* exclusive "or" operator node */
	    op = op_xor;
	binary:				/* common for all binary nodes */
	    describe_tree(tree->tr_b.tb_left, &left);
	    describe_tree(tree->tr_b.tb_right, &right);
	    bu_vls_putc(str, '(');
	    bu_vls_vlscatzap(str, &left);
	    bu_vls_printf(str, " %c ", op);
	    bu_vls_vlscatzap(str, &right);
	    bu_vls_putc(str, ')');
	    break;
	case OP_NOT:
	    bu_vls_strcat(str, "(!");
	    describe_tree(tree->tr_b.tb_left, str);
	    bu_vls_putc(str, ')');
	    break;
	case OP_GUARD:
	    bu_vls_strcat(str, "(G");
	    describe_tree(tree->tr_b.tb_left, str);
	    bu_vls_putc(str, ')');
	    break;
	case OP_XNOP:
	    bu_vls_strcat(str, "(X");
	    describe_tree(tree->tr_b.tb_left, str);
	    bu_vls_putc(str, ')');
	    break;
	case OP_NOP:
	    bu_vls_strcat(str, "NOP");
	    break;
	default:
	    bu_exit(1, "ERROR: describe_tree() got unrecognized op (%d)\n", tree->tr_op);
    }
}



/**
 * @brief This routine is called when a region is first encountered in the
 * hierarchy when processing a tree
 *
 *      @param tsp tree state (for parsing the tree)
 *      @param pathp A listing of all the nodes traversed to get to this node in the database
 *      @param combp the combination record for this region
 */
int
region_start(struct db_tree_state *tsp,
	     const struct db_full_path *pathp,
	     const struct rt_comb_internal *combp,
	     void *client_data)
{
    char *name;
    struct directory *dp;
    struct bu_vls str = BU_VLS_INIT_ZERO;
    struct user_data *your_stuff = (struct user_data *)client_data;

    RT_CK_DBTS(tsp);

    name = db_path_to_string(pathp);
    bu_log("region_start %s\n", name);
    bu_free(name, "reg_start name");

    bu_log("data = %ld\n", your_stuff->data);
    rt_pr_tol(&your_stuff->tol);

    dp = DB_FULL_PATH_CUR_DIR(pathp);

    /* here is where the conversion should be done */
    if (combp->region_flag)
	printf("Write this region (name=%s) as a part in your format:\n", dp->d_namep);
    else
	printf("Write this combination (name=%s) as an assembly in your format:\n", dp->d_namep);

    describe_tree(combp->tree, &str);

    printf("\t%s\n\n", bu_vls_addr(&str));

    bu_vls_free(&str);

    return 0;
}


/**
 * @brief This is called when all sub-elements of a region have been processed by leaf_func.
 *
 *      @param tsp
 *      @param pathp
 *      @param curtree
 *
 *      @return TREE_NULL if data in curtree was "stolen", otherwise db_walk_tree will
 *      clean up the data in the union tree * that is returned
 *
 * If it wants to retain the data in curtree it can by returning TREE_NULL.  Otherwise
 * db_walk_tree will clean up the data in the union tree * that is returned.
 *
 */
union tree *
region_end (struct db_tree_state *tsp,
	    const struct db_full_path *pathp,
	    union tree *curtree,
	    void *UNUSED(client_data))
{
    char *name;

    RT_CK_DBTS(tsp);

    name = db_path_to_string(pathp);
    bu_log("region_end   %s\n", name);
    bu_free(name, "region_end name");

    return curtree;
}

/*declare torusflag*/
int flag = 0 ;
double Vadd[3];
double maga;
double magb;
double magc;
double magd;


/* This routine is called by the tree walker (db_walk_tree)
 * for every primitive encountered in the trees specified on the command line */
union tree *
primitive_func(struct db_tree_state *tsp,
	       const struct db_full_path *pathp,
	       struct rt_db_internal *ip,
	       void *UNUSED(client_data))
{
    int i;
    struct directory *dp;
    char *name;
    dp = DB_FULL_PATH_CUR_DIR(pathp);

    RT_CK_DBTS(tsp);

    name = db_path_to_string(pathp);
    bu_log("leaf_func    %s\n", name);
    bu_free(name, "region_end name");

    /* handle each type of primitive (see h/rtgeom.h) */
    if (ip->idb_major_type == DB5_MAJORTYPE_BRLCAD) {
	switch (ip->idb_type) {
	    /* most commonly used primitives */
	  
	    case ID_TOR:	/* torus */
		{ 
		    struct rt_tor_internal *tor = (struct rt_tor_internal *)ip->idb_ptr;
		    if ( flag == 0 )
		    {
		    torusmacro();
		    flag = 1;
		    }
		    printf(" \nobject {\tTorus (\n");
		    printf("\t< %g, %g, %g>, ", V3ARGS(tor->v));
		    printf("<%g, %g, %g>, ", V3ARGS(tor->h));
		    printf(" %g , ", tor->r_a);
		    printf("%g )", tor->r_h);
		    printf(" texture{ pigment{ LightBlue} }}\n");
		    break;
		}
	    case ID_TGC: /* truncated general cone frustum */
		{
		    /* this primitive includes circular cross-section
		     * cones and cylinders
		     */
	    struct rt_tgc_internal *tgc = (struct rt_tgc_internal *)ip->idb_ptr;
	    VADD2(Vadd, tgc->v,tgc->h);
	    maga = MAGNITUDE(tgc->a);
	    magb = MAGNITUDE(tgc->b);
	    magc = MAGNITUDE(tgc->c);
	    magd = MAGNITUDE(tgc->d);
	    if(EQUAL(MAGNITUDE(tgc->a), MAGNITUDE(tgc->c)))     /* Cylinder */
	    { 
		    printf("\tcylinder\n\t    {\n ");
		   	printf("\t<%g %g %g>,\n", V3ARGS(tgc->v));
		  	printf("\t<%g %g %g>,  ", Vadd[0], Vadd[1],Vadd[2]);
		    printf("%g\n", maga);
		    printf("\t    texture{ pigment{ lightblue } }}\n");
		}  
		else if(EQUAL(MAGNITUDE(tgc->a), MAGNITUDE(tgc->b))) /* Cone */
		{
			printf("\tCone\n\t    {\n ");
			printf("\t<%g %g %g>,  ", V3ARGS(tgc->v));
		  	printf("%g,\n", maga);
			printf("\t    <%g %g %g>,  ", V3ARGS(tgc->h));
			printf("%g\n", magc);
			printf("\t    texture{ pigment{ lightblue } }}\n");
		}  
		else
		{	printf("#include \"shapes.inc\"\n");
			printf("\tobject{ Supercone(\n");
			printf("\t<%g, %g, %g>,  ", V3ARGS(tgc->v));
		  	printf("%g, %g ,\n", maga, magb);
			printf("<%g, %g, %g>,  ", Vadd[0], Vadd[1],Vadd[2]);
			printf("%g, %g)", magc, magd);
			printf("\t    texture{ pigment{ color rgb<0.65,1,0> } }}\n");
		}  
		break;
		}
	    case ID_REC: /* right elliptical cylinder */
		{
		    /* This primitive includes circular cross-section
		     * cones and cylinders
		     */
		    struct rt_tgc_internal *tgc = (struct rt_tgc_internal *)ip->idb_ptr;
		    printf("Write this TGC (name=%s) in your format:\n", dp->d_namep);
		    printf("\tV=(%g %g %g)\n", V3ARGS(tgc->v));
		    printf("\tH=(%g %g %g)\n", V3ARGS(tgc->h));
		    printf("\tA=(%g %g %g)\n", V3ARGS(tgc->a));
		    printf("\tB=(%g %g %g)\n", V3ARGS(tgc->b));
		    printf("\tC=(%g %g %g)\n", V3ARGS(tgc->c));
		    printf("\tD=(%g %g %g)\n", V3ARGS(tgc->d));
		    break;
		}
        case ID_ELL:
        {
		    /* ellipsoids */
		    struct rt_ell_internal *ell = (struct rt_ell_internal *)ip->idb_ptr;
		    maga = MAGNITUDE(ell->a);
		    magb = MAGNITUDE(ell->b);
		    magc = MAGNITUDE(ell->c);
		    printf("#include \"shapes.inc\"\nobject{\n\t\tSpheroid(\n");
		    printf("\t<%g, %g, %g>,\n", V3ARGS(ell->v));
		    printf("< %g ,", magb);
		    printf(" %g ,", maga);
		    printf(" %g > )", magc);
		    printf(" pigment{ LightBlue}\n\t}\n");
		    break;
		}
	    case ID_SPH:
		{
		    /* spheres*/
		    struct rt_ell_internal *ell = (struct rt_ell_internal *)ip->idb_ptr;
		    printf("sphere{\n");
		    printf("\t<%g, %g, %g>,\n", V3ARGS(ell->v));
		    printf("\t %g \n//%g%g\n", V3ARGS(ell->a));
		    printf(" pigment{ LightBlue}\n\t}\n");
		    break;
		}
        case ID_ARB8:       /* convex primitive with from four to six faces */
		{
		    /* this primitive may have degenerate faces
		    * faces are: 0123, 7654, 0347, 1562, 0451, 3267
		    * (points listed above in counter-clockwise order)
		    */
		    struct rt_arb_internal *arb = (struct rt_arb_internal *)ip->idb_ptr;

		    char coordinates[] = {'b','c','h','g','a','d','e','f'};
		    printf("#declare Default_texture = pigment {rgb 0.8}\n");
		    for(i=0; i<8; i++)
		    {
		    printf("#declare %c = <%g, %g, %g>;\n", coordinates[i], V3ARGS(arb->pt[i]));
		    }
		    printf("#declare Box = mesh{\n");
		    printf("triangle{a,b,c}\ntriangle{a,c,d}\ntriangle{a,d,f}\n");
		    printf("triangle{e,d,f}\ntriangle{c,d,e}\ntriangle{c,e,h}\n");
		    printf("triangle{a,b,g}\ntriangle{a,f,g}\ntriangle{b,c,g}\n");
		    printf("triangle{g,h,c}\ntriangle{e,f,g}\ntriangle{e,g,h}\n");
		    printf("texture{Default_texture}\n}\n Box\n");
			break;
		}


		/* other primitives, left as an exercise to the reader */

		case ID_BOT:        /* Bag O' Triangles */
		case ID_ARS:
		    /* series of curves
		    * each with the same number of points
		    */
		case ID_HALF:
                     /* half universe defined by a plane */
		case ID_POLY:
		    /* polygons (up to 5 vertices per) */
		case ID_BSPLINE:
		   /* NURB surfaces */
		case ID_NMG:
		   /* N-manifold geometry */
		case ID_ARBN:
		case ID_DSP:
		   /* Displacement map (terrain primitive) */
		   /* the DSP primitive may reference an external file or binunif object */
		case ID_HF:
		   /* height field (terrain primitive) */
		   /* the HF primitive references an external file */
		case ID_EBM:
		   /* extruded bit-map */
		   /* the EBM primitive references an external file */
		case ID_VOL:
		   /* the VOL primitive references an external file */
		case ID_PIPE:
		case ID_PARTICLE:
		{
		    struct rt_part_internal *part = (struct rt_part_internal *)ip->idb_ptr;
		    printf("#include \"shapes.inc\"\n");
		    printf("object{\n\t Round_Cone2(\n");
		    printf("\t\t<%g %g %g>,", V3ARGS(part->part_V) );
		    printf(" %g,\n", part->part_vrad );
		    printf("\t\t <%g %g %g>,", V3ARGS(part->part_H) );
		    printf(" %g, 0)\n}\n", part->part_hrad );
		    break;
		}
		case ID_RPC:
		case ID_RHC:
		case ID_EPA:
		case ID_EHY:
		case ID_ETO:
		case ID_GRIP:
		case ID_SKETCH:
		case ID_EXTRUDE:
		    /* note that an extrusion references a sketch, make sure you convert
		    * the sketch also
		    */
	    

		/* other primitives, left as an exercise to the reader */

	    
	    default:
		bu_log("Primitive %s is an unsupported or unrecognized type (%d)\n", dp->d_namep, ip->idb_type);
		break;
	}
    } else {
	switch (ip->idb_major_type) {
	    case DB5_MAJORTYPE_BINARY_UNIF:
		{
		    /* not actually a primitive, just a block of storage for data
		     * a uniform array of chars, ints, floats, doubles, ...
		     */
		    struct rt_binunif_internal *bin = (struct rt_binunif_internal *)ip->idb_ptr;

		    if (bin)
			printf("Found a binary object (%s)\n\n", dp->d_namep);
		    break;
		}
	    default:
		bu_log("Major type of %s is unrecognized type (%d)\n", dp->d_namep, ip->idb_major_type);
		break;
	}
    }

    return (union tree *) NULL;
}


int
main(int argc, char *argv[])
{
    static const char usage[] = "Usage: %s [-xX lvl] [-a abs_tol] [-r rel_tol] [-n norm_tol] [-o out_file] [-C Camera_loc] [-V Look_at] [-L Light_loc] [-l Light_col] [-D default] brlcad_db.g object(s)\n";

    struct user_data your_data = {0, BN_TOL_INIT_ZERO};
    int i;
    int c;
    char idbuf[132] = {0};

    struct rt_i *rtip;
    struct db_tree_state init_state;

    bu_setprogname(argv[0]);
    bu_setlinebuf(stderr);

    /* calculational tolerances
     * mostly used by NMG routines
     */
    your_data.tol.magic = BN_TOL_MAGIC;
    your_data.tol.dist = 0.0005;
    your_data.tol.dist_sq = your_data.tol.dist * your_data.tol.dist;
    your_data.tol.perp = 1e-6;
    your_data.tol.para = 1 - your_data.tol.perp;
    

    /* Get command line arguments. */
    while ((c = bu_getopt(argc, argv, "t:a:n:o:r:x:X:C:V:L:l:c:D")) != -1) {
	float a1, a2, a3, a4, b1, b2, b3, b4,  c2, c3, c4;
	switch (c) {
	    case 't':		/* calculational tolerance */
		your_data.tol.dist = atof(bu_optarg);
		your_data.tol.dist_sq = your_data.tol.dist * your_data.tol.dist;
	    case 'o':		/* Output file name */
		/* grab output file name */
		break;
	    case 'x':		/* librt debug flag */
		sscanf(bu_optarg, "%x", &RTG.debug);
		bu_printb("librt RT_G_DEBUG", RT_G_DEBUG, DEBUG_FORMAT);
		bu_log("\n");
		break;
	    case 'X':		/* NMG debug flag */
		sscanf(bu_optarg, "%x", &RTG.NMG_debug);
		bu_printb("librt RTG.NMG_debug", RTG.NMG_debug, NMG_DEBUG_FORMAT);
		bu_log("\n");
		break;
	    case 'C':
		if(sscanf(bu_optarg, "%g, %g", &a1, &b1))
		{
			printf("the value of a b c %g, %g ", a1, b1);
			bu_log("\n");
		}
		else
		{
		        bu_exit(1, usage, argv[0]);
		}
		break;
	    case 'V':
		sscanf(bu_optarg, "%g %g %g", &a2, &b2, &c2);
		printf("Camera View point ");
		bu_log("\n");
		break;
	    case 'L':
		sscanf(bu_optarg, "%g %g %g", &a3, &b3, &c3);
		printf("Light");
		bu_log("\n");
		break;
	    case 'l':
		sscanf(bu_optarg, "%g %g %g", &a4, &b4, &c4);
		printf("Light colour");
		bu_log("\n");
		break;
	    case 'D':
		
		printf("\n#include\"colors.inc\"\n");
		printf("\nbackground { color Black }\n");
		printf("camera\n\t{\n\t\tlocation <0, 0, 40>\n\t\tlook_at <0, 0, 0>\n\t\t\t}\n");
		printf("light_source\n\t{\n\t\t<0, 0, 40> White\n\t\t}\n");
		bu_log("\n");
		break;
	    default:
		bu_exit(1, usage, argv[0]);
		break;
	}
    }

    if (bu_optind+1 >= argc) {
	bu_exit(1, usage, argv[0]);
    }

    /* Open BRL-CAD database */
    /* Scan all the records in the database and build a directory */
    rtip=rt_dirbuild(argv[bu_optind], idbuf, sizeof(idbuf));
    if (rtip == RTI_NULL) {
	bu_exit(1, "g-xxx: rt_dirbuild failure\n");
    }

    init_state = rt_initial_tree_state;

    bu_optind++;

    /* Walk the trees named on the command line
     * outputting combinations and primitives
     */
    for (i=bu_optind; i<argc; i++) {

	
	db_walk_tree(rtip->rti_dbip, argc - i, (const char **)&argv[i], 1 /* bu_avail_cpus() */,
		     &init_state, region_start, region_end, primitive_func, (void *) &your_data);
    }

    return 0;
}


void torusmacro()
{


		   { 
		    printf("#include\"transforms.inc\"\n");
		    printf("#macro Torus(Center, Normal, Radius1, Radius2)\n");
		    printf("\t torus{ Radius1, Radius2 Reorient_Trans(y, Normal) translate Center }\n#end\n\n");
		    }

}


/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
