use strict;

# input: kt_self.txt, matrix.txt
# output: plotly HTML
my $input='tax';
my $H=50;
my $n=0;
my @x; # labels
my @y; # labels
my @z; # csv row-vectors
if($input eq 'kt') {
  open(IN,"<kt_self.txt");
  while(<IN>) {
	next if(/^#/);
	$n++;
	last if($n>$H);
	chomp;
	my(@tmp)=split(/\t/);
        my($species)=pop(@tmp);
	push(@x,$species);
  }
} else { # tax.tab
  open(IN,"<tax.tab");
  while(<IN>) {
        next if(/^#/);
        my(@tmp)=split(/\t/);
        my($species)=$tmp[6];
	next if($species eq 'species');
        $n++;
        last if($n>$H);
        push(@x,$species);
  }
}
close(IN);
open(IN,"<matrix.txt");
while(<IN>) {
	next if(/qpid/);
	chomp;
	my($qpid,$csv)=split(/\t/);
	push(@y,$qpid);
	push(@z,$csv);
}
close(IN);
warn "# data sizes: x $#x y $#y z $#z\n";
# output HTML 
my $zstring=join("],\n[",@z);
my $ystring=join("\",\n\"",@y);
my $xstring=join("\",\n\"",@x);
my $height=5*($#y+1); 
if($height>32000) { $height=32000; }
if($height<1000) { $height=1000; }
print<<EOB;
<head>
  <!-- Plotly.js -->
  <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
</head>

<body>
Browser compatibility: Chrome recommended. Windows Edge does not display Plotly heatmap correctly.  
  <div id="myDiv"><!-- Plotly chart will be drawn inside this DIV --></div>
  <script>
    <!-- JAVASCRIPT CODE GOES HERE -->
var layout = {
  autosize: false,
  width: 1500,
  height: $height,
  margin: {
    l: 160,
    r: 100,
    b: 400,
    t: 100,
    pad: 4
  },
}
var data = [
  {
z: [ [$zstring] ],
y: [ "$ystring" ],
x: [ "$xstring" ],
colorscale: 'Picnic',
type: 'heatmap'
  }
];

Plotly.plot('myDiv', data, layout);

  </script>
</body>
EOB

