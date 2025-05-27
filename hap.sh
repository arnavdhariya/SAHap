#!/bin/bash
################## SKELETON: DO NOT TOUCH THESE 2 LINES
EXEDIR=`dirname "$0"`; BASENAME=`basename "$0" .sh`; TAB='	'; NL='
'
#################### ADD YOUR USAGE MESSAGE HERE, and the rest of your code after END OF SKELETON ##################
USAGE="USAGE: $BASENAME bla bla bla
PURPOSE: describe purpose in words"

################## SKELETON: DO NOT TOUCH CODE HERE
# check that you really did add a usage message above
USAGE=${USAGE:?"$0 should have a USAGE message before sourcing skel.sh"}
die(){ echo "$USAGE${NL}FATAL ERROR in $BASENAME:" "$@" >&2; exit 1; }
[ "$BASENAME" == skel ] && die "$0 is a skeleton Bourne Shell script; your scripts should source it, not run it"
echo "$BASENAME" | grep "[ $TAB]" && die "Shell script names really REALLY shouldn't contain spaces or tabs"
[ $BASENAME == "$BASENAME" ] || die "something weird with filename in '$BASENAME'"
warn(){ (echo "WARNING: $@")>&2; }
not(){ if eval "$@"; then return 1; else return 0; fi; }
newlines(){ awk '{for(i=1; i<=NF;i++)print $i}' "$@"; }
parse(){ awk "BEGIN{print $*}" </dev/null; }
which(){ echo "$PATH" | tr : "$NL" | awk '!seen[$0]{print}{++seen[$0]}' | while read d; do eval /bin/ls $d/$N; done 2>/dev/null | newlines; }

unset TMPDIR

[ "$MYTMP" ] || export MYTMP="/tmp"
export TMPDIR=${TMPDIR:-`mktemp -d $MYTMP/$BASENAME.XXXXXX`}
 trap "/bin/rm -rf $TMPDIR; exit" 0 1 2 3 15 # call trap "" N to remove the trap for signal N

#################### END OF SKELETON, ADD YOUR CODE BELOW THIS LINE

# a WIF file is one read per line
sed 's/ : #.*//' "$1" | # remove the final colon and comment text
    awk '{for(c=1;c<NF;c+=5)print $c,$(c+1),FNR}' | sort | # site, letter, and read, sorted by site
    gawk '{
	if(prevSL==$1" "$2) printf " %d", $3; # print read number of the site if site+letter agree with previous line
	else {
	    sameSL=0; if(start++) print ""; # print newline except first time
	    prevSL=$1" "$2; printf "%s %s\t%s",$1,$2,$3 # otherwise print new site, letter, and read
	}
    }' | tee $TMPDIR/slRs.txt | # at this point we have "site letter {set of reads that have this letter at this site}"
    hawk '{for(i=3;i<NF;i++) for(j=i+1;j<=NF;j++) ++agree[MIN($i,$j)][MAX($i,$j)]} # accumulate #sites at which reads agree
	END{for(i in agree)for(j in agree[i]) print agree[i][j],i,j}' |
    sort -nr | tee $TMPDIR/nrr.txt | # countAgree r1 r2
    hawk '{
	    r1=$2;r2=$3;
	    if(!(r1 in group)&&!(r2 in group)) { # neither are in a group
		group[r1]=group[r2]=++numGroups; set[numGroups][r1]=set[numGroups][r2]=$1;
	    }
	    else if((r1 in group)&&(r2 in group)) { # both are in a group
		if(group[r1] != group[r2]) { # merge the groups... merge them into the bigger scoring one
		    ASSERT(set[group[r1]][r1]>0 && set[group[r2]][r2]>0);
		    if(set[group[r1]][r1] > set[group[r2]][r2]) {src=r2;dest=r1}
		    else                                        {src=r1;dest=r2}
		    for(r in set[group[src]]) {set[group[dest]][r]=set[group[dest]][dest]; group[r]=group[dest]}
		    delete set[group[src]];
		}
	    } # at this point EXACTLY one is in a group
	    else if(r1 in group){group[r2]=group[r1]; set[group[r1]][r2]=$1}
	    else if(r2 in group){group[r1]=group[r2]; set[group[r2]][r1]=$1}
	}
	END{for(g in set){printf "%s\t", g; for(r in set[g])printf " %s",r; print ""}}'
