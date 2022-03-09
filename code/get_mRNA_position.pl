#北京组学生物科技有限公司
#学习perl语言：
#perl入门：https://study.163.com/course/introduction/1006448023.htm?share=1&shareId=1030291076
#perl高级：https://study.163.com/course/introduction/1004833023.htm?share=1&shareId=1030291076

die "perl $0 <mRNAid.list>  <gff3>  <out.gff>" unless(@ARGV==3);
use Data::Dumper;


open IN,"$ARGV[0]" or die "$!; can't open file $ARGV[0]\n";

open OUT,">$ARGV[2]" or die "$!; can't open file $ARGV[2]\n";

my %mRNAID=();

while(<IN>){
	chomp;
	my@tmp=split(/\s+/);
	$mRNAID{$tmp[0]}=1;

}

#print Dumper(\%mRNAID);
close(IN);

open IN,"$ARGV[1]" or die "$!; can't open file $ARGV[1]\n";
while(<IN>){
	chomp;

	next if /^#/;
	@tmp=split(/\t/);
	#if($tmp[2]=~/gene/ && $tmp[0]=~/^\d+/ && $tmp[-1]=~/protein_coding/){
	if($tmp[2] =~/^mRNA$/i or $tmp[2] =~/^transcript$/i){
		my($id)=($tmp[-1]=~/ID=([^;]+)/);
		if (exists $mRNAID{$id}) {
			print OUT "$tmp[0]\t$id\t$tmp[3]\t$tmp[4]\n";
		}
		
	}
}

close(IN);
close(OUT);
