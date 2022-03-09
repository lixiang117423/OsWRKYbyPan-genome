#北京组学生物科技有限公司
#学习perl语言：
#perl入门：https://study.163.com/course/introduction/1006448023.htm?share=1&shareId=1030291076
#perl高级：https://study.163.com/course/introduction/1004833023.htm?share=1&shareId=1030291076


die "perl $0 <in>  <out>" unless(@ARGV==2);
use Bio::SeqIO;
use Bio::Seq;
use Bio::Tools::SeqStats;
use Bio::Tools::pICalculator;
use Data::Dumper;
#读入序列
my $in = Bio::SeqIO->new(
	-file   => "$ARGV[0]",
	-format => 'Fasta'
);

open OUT,">$ARGV[1]" or die "$!; can't open file $ARGV[1]\n";
print OUT "#ID\tlength\tMW(Da)\tpI\n";
my $calc = Bio::Tools::pICalculator->new(-places => 2,-pKset => 'EMBOSS');


#逐条读取序列并计算
while ( my $seq = $in->next_seq() ) {
	#my ( $id, $sequence, $desc ) = ( $seq->id, $seq->seq, $seq->desc );
	my $weight = Bio::Tools::SeqStats ->get_mol_wt($seq);
	$calc->seq($seq);
    my $iep = $calc->iep;
    print OUT sprintf("%s\t%s\t%s\t%s\n",
                  $seq->id,
                  $seq->length,
                  "$weight->[0]",
                  $iep);
}
$in->close();
close(OUT);
