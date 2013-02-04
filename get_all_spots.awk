BEGIN{
	printf "[ "
	if(s == "")
		skip=1
	else
		skip = s+0
}

FNR==1 { 
	ignore_next=skip
	print FILENAME > "/dev/stderr"
}
/MT19/{next}
/INTERMEDIATE/{next}

/MAIN/{
	$1=""
}

/PASS0/{next}
/PASS./{$1 = ""}

NF>2{
	if(ignore_next<=0)
		printf " "$0" "
	ignore_next--
}

END{
	print " ]"
}
