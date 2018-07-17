#!/usr/bin/env perl

# Copyright (c) 2018 Hikoyu Suzuki
# This software is released under the MIT License.

use strict;
use warnings;
use Getopt::Std;
#use threads;

# ソフトウェアを定義
### 編集範囲 開始 ###
my $software = "teana.pl";	# ソフトウェアの名前
my $version = "ver.1.0.1";	# ソフトウェアのバージョン
my $note = "TEANA is Terminus Extending Assmbler with Nearby Alignment.\n  This software assembles short reads aligned to the both ends of seed sequences and extends them.";	# ソフトウェアの説明
my $usage = "<seed.fa> <in1_r1.fq[,in2_r1.fq,...]> <in1_r2.fq[,in2_r2.fq,...]>";	# ソフトウェアの使用法 (コマンド非使用ソフトウェアの時に有効)
### 編集範囲 終了 ###

# コマンドを定義
my %command;
### 編集範囲 開始 ###
# コマンドを追加
### 編集範囲 終了 ###
my @command_list = sort(keys(%command));

# 指定されたコマンドを確認
my $specified_command = shift(@ARGV) if @command_list and @ARGV;
&exception::error("unknown command: $specified_command") if $specified_command and !grep {$_ eq $specified_command} @command_list;

# 共通オプションを定義
my %option;
### 編集範囲 開始 ###
$option{"v"} = "\tKeep intermediate files for assembly";
$option{"o PATH "} = "Path to output directory [.]";
$option{"p STR "} = "Presets for aligning reads <very-fast|fast|sensitive|very-sensitive> [sensitive]";
$option{"c INT "} = "Maximum clipping size <0-> [0]";
$option{"d INT "} = "Minimum k-mer depth <1-> [1]";
$option{"j INT "} = "Minimum abundance of junction reads <1-> [2]";
$option{"k INT "} = "K-mer size for assembly <19-31> [25]";
$option{"l INT "} = "Minimum fragment length <1-> [200]";
$option{"m INT "} = "Maximum fragment length <1-> [800]";
$option{"n INT "} = "Cutoff overlap length <0-> [300]";
$option{"t INT "} = "Number of parallel worker threads <1-> [1]";
$option{"i FLOAT "} = "Cutoff overlap identity <0-1> [0.99]";
# オプションを追加
### 編集範囲 終了 ###

# コマンドごとのオプション定義を取得
&{\&{"${specified_command}::define"}} if $specified_command;
my @option_list = sort(keys(%option));

# ヘルプを表示 (引数未指定時)
&exception::help if !@ARGV and !-p STDIN;

# オプションの入力処理
my %opt;
$_ = join("", @option_list);
$_ =~ s/\s+\S+\s+/:/g;
getopts($_, \%opt);

# 未指定オプションのデフォルト値を入力
foreach (@option_list) {
	$opt{substr($_, 0, 1)} = substr($option{$_}, index($option{$_}, "[") + 1, index($option{$_}, "]") - index($option{$_}, "[") - 1) if $option{$_} =~ /\[.+\]$/ and !defined($opt{substr($_, 0, 1)});
}

### 編集範囲 開始 ###
# 追加のモジュールを宣言
use List::Util;

# 処理を追加
### 編集範囲 終了 ###

# メインルーチンを実行
&main;
exit(0);

# メインルーチン
sub main {
	### 編集範囲 開始 ###
	# 指定された共通オプションを確認
	&exception::error("unknown presets specified: -p $opt{p}") if $opt{"p"} ne "very-fast" and $opt{"p"} ne "fast" and $opt{"p"} ne "sensitive" and $opt{"p"} ne "very-sensitive";
	&exception::error("specify INT >= 0: -c $opt{c}") if $opt{"c"} !~ /^\d+$/;
	&exception::error("specify INT >= 1: -d $opt{d}") if $opt{"d"} !~ /^\d+$/ or $opt{"d"} < 1;
	&exception::error("specify INT >= 1: -j $opt{j}") if $opt{"j"} !~ /^\d+$/ or $opt{"j"} < 1;
	&exception::error("specify 19 <= INT <= 32: -k $opt{k}") if $opt{"k"} !~ /^\d+$/ or $opt{"k"} < 19 or $opt{"k"} > 32;
	&exception::error("specify INT >= 1: -l $opt{l}") if $opt{"l"} !~ /^\d+$/ or $opt{"l"} < 1;
	&exception::error("specify INT >= $opt{l}: -m $opt{m}") if $opt{"m"} !~ /^\d+$/ or $opt{"m"} < $opt{"l"};
	&exception::error("specify INT >= 0: -n $opt{n}") if $opt{"n"} !~ /^\d+$/;
	&exception::error("specify INT >= 1: -t $opt{t}") if $opt{"t"} !~ /^\d+$/ or $opt{"t"} < 1;
	&exception::error("specify FLOAT 0-1: -i $opt{i}") if $opt{"i"} !~ /^\d+$|^\d+\.\d+$|^\d+[eE]-?\d+$|^\d+\.\d+[eE]-?\d+$/ or $opt{"i"} > 1;
	
	# 依存関係を確認
	my $software_list = &check_dependencies(["bowtie2", "Bridger.pl", "makeblastdb", "blastn", "blastdbcmd"]);
	&exception::error("bowtie2 not installed") if !$software_list->{"bowtie2"};
	&exception::error("bridger not installed") if !$software_list->{"Bridger.pl"};
	&exception::error("makeblastdb not installed") if !$software_list->{"makeblastdb"};
	&exception::error("blastn not installed") if !$software_list->{"blastn"};
	&exception::error("blastdbcmd not installed") if !$software_list->{"blastdbcmd"};
	
	# 入力ファイルを取得
	my ($seed_file, $read1_file, $read2_file) = @ARGV;
	
	# 入力ファイルを確認
	&exception::error("seed file not specified") if !defined($seed_file);
	&exception::error("read1 files not specified") if !defined($read1_file);
	&exception::error("read2 files not specified") if !defined($read2_file);
	&check_files([$seed_file, split(/,/, "$read1_file,$read2_file")]);
	
	# 出力先ディレクトリを作成
	mkdir($opt{"o"}) or &exception::error("failed to make directory: $opt{o}") if !-d $opt{"o"};
	
	# 処理を定義
	my $build_index = "bowtie2-build --threads $opt{t} $seed_file $opt{o}/seed >$opt{o}/bowtie2_seed.log 2>&1";
	my $align_reads = "bowtie2 --end-to-end --fr -p $opt{t} -I $opt{l} -X $opt{m} --$opt{p} -x $opt{o}/seed -1 $read1_file -2 $read2_file 2>>$opt{o}/bowtie2_seed.log";
	my $assemble_reads = "Bridger.pl --CPU $opt{t} --seqType fa --SS_lib_type FR --min_kmer_coverage $opt{d} --min_junction_coverage $opt{j} --kmer_length $opt{k}";
	my $make_database = "makeblastdb -in $seed_file -out $opt{o}/seed -dbtype nucl -parse_seqids";
	my $search_sequence = "blastn -num_threads $opt{t} -task megablast -strand plus -db $opt{o}/seed -outfmt '6 qseqid qstart qend qlen sstart send slen nident length'";
	my $obtain_seed = "blastdbcmd -db $opt{o}/seed -dbtype nucl -outfmt %s";
	
	# シード配列のbowtie2インデックスを作成
	print STDERR "Creating an index of seeds...";
	system($build_index) and &exception::error("failed to create an index: $seed_file");
	print STDERR "completed\n";
	
	# 変数を宣言
	my %seed_len = ();
	my %read_files = ();
	my @buf = ("");
	my $min_read_len = 0;
	my $num_seeds = 0;
	
	# シードサマリーファイルを作成
	open(SEED_SUMMARY, ">", "$opt{o}/seed_summary.txt") or &exception::error("failed to make file: $opt{o}/seed_summary.txt");
	
	# アラインメントを実行
	open(ALIGNMENT_OUT, "-|", $align_reads) or &exception::error("failed to align reads: $read1_file & $read2_file -> $seed_file");
	
	# アラインメントの出力を読み込みながら処理
	print STDERR "Aligning reads to seeds...";
	while (my $line = <ALIGNMENT_OUT>) {
		# タブ文字で分割
		my @col = split(/\t/, $line);
		
		# シード配列長を取得
		$seed_len{substr($col[1], 3)} = substr($col[2], 3) if $col[0] eq '@SQ';
		
		# ヘッダー行の場合は除く
		next if $line =~ /^@/;
		
		# フラグ条件を満たす場合は除く
		next if $col[1] & 0xF02;
		
		# 逆鎖にマップされている場合はリード配列を相補鎖に変換
		&complementary($col[9]) if $col[1] & 0x010;
		
		# データバッファーに相方の情報が存在する場合
		if ($col[0] eq $buf[0] and (($col[1] | $buf[1]) & 0x0C0) == 0x0C0) {
			# 自身がシード配列の上流側の逆鎖にマップされている場合
			if (($col[1] & 0x014) == 0x010 and $col[3] <= $opt{"m"}) {
				# リード配列ファイルを作成
				&create_read_file($read_files{$col[2]}, $num_seeds) and print SEED_SUMMARY "seed$num_seeds\t$col[2]\n" and system("echo '$col[2]' >$opt{o}/seed$num_seeds.txt 2>/dev/null") and &exception::error("failed to make file: seed$num_seeds.txt");
				
				# リード配列ファイルにfasta形式で出力
				print {$read_files{$col[2]}->[0]} ">$buf[0]\n$buf[9]\n";
				print {$read_files{$col[2]}->[1]} ">$col[0]\n$col[9]\n";
			}
			
			# 自身がシード配列の下流側の順鎖にマップされている場合
			elsif (($col[1] & 0x014) == 0x000 and $col[3] + length($col[9]) > $seed_len{$col[2]} - $opt{"m"}) {
				# リード配列ファイルを作成
				&create_read_file($read_files{$col[2]}, $num_seeds) and print SEED_SUMMARY "seed$num_seeds\t$col[2]\n" and system("echo '$col[2]' >$opt{o}/seed$num_seeds.txt 2>/dev/null") and &exception::error("failed to make file: seed$num_seeds.txt");
				
				# リード配列ファイルにfasta形式で出力
				print {$read_files{$col[2]}->[2]} ">$col[0]\n$col[9]\n";
				print {$read_files{$col[2]}->[3]} ">$buf[0]\n$buf[9]\n";
			}
			
			# 相方がシード配列の上流側の逆鎖にマップされている場合
			if (($buf[1] & 0x014) == 0x010 and $buf[3] <= $opt{"m"}) {
				# リード配列ファイルを作成
				&create_read_file($read_files{$buf[2]}, $num_seeds) and print SEED_SUMMARY "seed$num_seeds\t$buf[2]\n" and system("echo '$buf[2]' >$opt{o}/seed$num_seeds.txt 2>/dev/null") and &exception::error("failed to make file: seed$num_seeds.txt");
				
				# リード配列ファイルにfasta形式で出力
				print {$read_files{$buf[2]}->[0]} ">$col[0]\n$col[9]\n";
				print {$read_files{$buf[2]}->[1]} ">$buf[0]\n$buf[9]\n";
			}
			
			# 相方がシード配列の下流側の順鎖にマップされている場合
			elsif (($buf[1] & 0x014) == 0x000 and $buf[3] + length($buf[9]) > $seed_len{$buf[2]} - $opt{"m"}) {
				# リード配列ファイルを作成
				&create_read_file($read_files{$buf[2]}, $num_seeds) and print SEED_SUMMARY "seed$num_seeds\t$buf[2]\n" and system("echo '$buf[2]' >$opt{o}/seed$num_seeds.txt 2>/dev/null") and &exception::error("failed to make file: seed$num_seeds.txt");
				
				# リード配列ファイルにfasta形式で出力
				print {$read_files{$buf[2]}->[2]} ">$buf[0]\n$buf[9]\n";
				print {$read_files{$buf[2]}->[3]} ">$col[0]\n$col[9]\n";
			}
			
			# 最小リード長を更新
			$min_read_len = length($col[9]) + length($buf[9]) if !$min_read_len or length($col[9]) + length($buf[9]) < $min_read_len;
		}
		
		# データバッファーを更新
		@buf = @col;
	}
	print STDERR "completed\n";
	
	# アラインメントの出力を閉じる
	close(ALIGNMENT_OUT);
	
	# シードサマリーファイルを閉じる
	close(SEED_SUMMARY);
	
	# リード配列ファイルを閉じる
	foreach my $file (values(%read_files)) {
		close($file->[0]);
		close($file->[1]);
		close($file->[2]);
		close($file->[3]);
	}
	
	# シード配列のBLAST+データベースを作成
	print STDERR "Creating a database of seeds...";
	system("$make_database >$opt{o}/makeblastdb_seed.log 2>&1") and &exception::error("failed to create a database: $seed_file");
	print STDERR "completed\n";
	
	# 変数を宣言
	my @region = ("upstream", "downstream");
	my $pair_distance = $opt{"m"} - $min_read_len;
	
	# 延長コンティグファイルを作成
	open(EXTENDED_CONTIGS, ">", "$opt{o}/extended_contigs.fa") or &exception::error("failed to make file: $opt{o}/extended_contigs.fa");
	
	# 各シード配列について処理
	foreach my $i (1..$num_seeds) {
		# 変数を宣言
		my @contig_seq = ({"" => ""}, {"" => ""});
		my @candidate_contigs = ({}, {});
		
		# シード配列を出力
		open(SEED, "-|", "$obtain_seed -entry_batch $opt{o}/seed$i.txt") or &exception::error("failed to obtain seed$i sequences");
		
		# シード配列を読み込む
		my $seed_seq = <SEED>;
		
		# 改行コードを除く
		chomp($seed_seq);
		
		# シード配列の出力を閉じる
		close(SEED);
		
		# 上流側と下流側それぞれについて処理
		foreach my $j (0..1) {
			# bridgerによるアセンブルを実行
			print STDERR "Assembling reads of seed$i $region[$j] region...";
			system("$assemble_reads --pair_gap_length $pair_distance --left $opt{o}/seed${i}_$region[$j]_r1.fa --right $opt{o}/seed${i}_$region[$j]_r2.fa --output $opt{o}/seed${i}_$region[$j] >$opt{o}/Bridger_seed${i}_$region[$j].log 2>&1") and &exception::caution("failed to assemble reads") or print STDERR "completed\n";
			
			# コンティグファイルを移動
			system("mv $opt{o}/seed${i}_$region[$j]/Bridger.fasta $opt{o}/seed${i}_$region[$j]_contigs.fa") and &exception::error("failed to move file: $opt{o}/seed${i}_$region[$j]/Bridger.fasta -> $opt{o}/seed${i}_$region[$j]_contigs.fa") if -f "$opt{o}/seed${i}_$region[$j]/Bridger.fasta";
			
			# 中間ファイルを削除 (-v未指定時)
			system("rm -rf $opt{o}/seed${i}_$region[$j]/") and &exception::error("failed to remove files: $opt{o}/seed${i}_$region[$j]/") if !$opt{"v"} and -d "$opt{o}/seed${i}_$region[$j]/";
			
			# アセンブルに失敗した場合は除く
			next if !-f "$opt{o}/seed${i}_$region[$j]_contigs.fa";
			
			# 変数を宣言
			my $contig_id = "";
			
			# コンティグファイルを開く
			open(CONTIG, "<", "$opt{o}/seed${i}_$region[$j]_contigs.fa") or &exception::error("failed to open file: $opt{o}/seed${i}_$region[$j]_contigs.fa");
			
			# コンティグファイルを読み込みながら処理
			print STDERR "Loading contigs of seed$i $region[$j] region...";
			while (my $line = <CONTIG>) {
				# 改行コードを除く
				chomp($line);
				
				# コンティグIDを取得
				($contig_id) = split(/\s+/, substr($line, 1)) and next if $line =~ /^>/;
				
				# コンティグ配列を取得
				$contig_seq[$j]->{$contig_id} .= $line;
			}
			print STDERR "completed\n";
			
			# コンティグファイルを閉じる
			close(CONTIG);
			
			# 相同性検索を実行
			open(HOMOLOGY_SEARCH_OUT, "-|", "$search_sequence -query $opt{o}/seed${i}_$region[$j]_contigs.fa -seqidlist $opt{o}/seed$i.txt");
			
			# 相同性検索の出力を読み込みながら処理
			print STDERR "Detecting overlaps in seed$i $region[$j] region...";
			while (<HOMOLOGY_SEARCH_OUT>) {
				# タブ文字で分割
				my @col = split(/\t/);
				
				# 連結条件を満たす場合は候補に追加
				(vec($candidate_contigs[$j]->{$col[0]}, 0, 32), vec($candidate_contigs[$j]->{$col[0]}, 1, 32)) = ($col[$j + 1] - 1, $col[$j + 4] - 1) if $col[2 + $j * 3] >= $col[3 + $j * 3] - $opt{"c"} and $col[4 - $j * 3] <= 1 + $opt{"c"} and $col[7] / $col[8] >= $opt{"i"} and $col[8] >= $opt{"n"};
			}
			print STDERR "completed\n";
			
			# 相同性検索の出力を閉じる
			close(HOMOLOGY_SEARCH_OUT);
			
			# 連結条件を満たすコンティグが得られなかった場合はデフォルト値を設定
			(vec($candidate_contigs[$j]->{""}, 0, 32), vec($candidate_contigs[$j]->{""}, 1, 32)) = (0, $j * length($seed_seq)) if !%{$candidate_contigs[$j]};
		}
		
		# アセンブルされた配列とシード配列を連結して出力
		print STDERR "Gluing seed$i...";
		foreach my $upstream (sort(keys(%{$candidate_contigs[0]}))) {
			my $tmp_seq = substr($contig_seq[0]->{$upstream}, 0, vec($candidate_contigs[0]->{$upstream}, 0, 32)) . substr($seed_seq, vec($candidate_contigs[0]->{$upstream}, 1, 32));
			foreach my $downstream (sort(keys(%{$candidate_contigs[1]}))) {
				print EXTENDED_CONTIGS ">$upstream:seed$i:$downstream\n";
				print EXTENDED_CONTIGS substr($tmp_seq, 0, vec($candidate_contigs[0]->{$upstream}, 0, 32) + vec($candidate_contigs[1]->{$downstream}, 1, 32) - vec($candidate_contigs[0]->{$upstream}, 1, 32)) . substr($contig_seq[1]->{$downstream}, vec($candidate_contigs[1]->{$downstream}, 0, 32)), "\n";
			}
		}
		print STDERR "completed\n";
	}
	
	# 延長コンティグファイルを閉じる
	close(EXTENDED_CONTIGS);
	### 編集範囲 終了 ###
	
	# コマンドの実行 (コマンド指定時)
	if ($specified_command) {&{\&{"${specified_command}::body"}};}
	
	### 編集範囲 開始 ###
	# 処理を追加
	### 編集範囲 終了 ###
	return(1);
}

### 編集範囲 開始 ###
# 相補鎖変換 complementary(配列)
sub complementary {
	# 引数を取得
	my ($seq) = @_;
	
	# 配列を逆順に並べ替える
	$seq = reverse($seq);
	
	# 相補的な塩基に置換
	$seq =~ tr/ATGCRYKMDBVH/TACGYRMKHVBD/;
	
	# 引数の配列を変更
	$_[0] = $seq;
	return(1);
}

# 依存関係確認 check_dependencies(コマンドリストリファレンス)
sub check_dependencies {
	# 引数を取得
	my ($command_list) = @_;
	
	# 変数を宣言
	my %result = ();
	
	# 各コマンドについて処理
	foreach my $command (@{$command_list}) {
		next if !$command;
		$result{$command} .= `$command -version 2>/dev/null`;
		$result{$command} .= `$command --version 2>/dev/null`;
	}
	
	# 結果を返す
	return(\%result);
}

# ファイル確認 check_files(ファイルリストリファレンス)
sub check_files {
	# 引数を取得
	my ($file_list) = @_;
	
	# 各ファイルについて処理
	foreach my $file (@{$file_list}) {
		next if !$file;
		&exception::error("file not found: $file") if !-f $file;
		&exception::error("file unreadable: $file") if !-r $file;
		&exception::error("null file specified: $file") if !-s $file;
	}
	return(1);
}

# リード配列ファイルを作成 create_read_file(ファイルハンドル, シードID)
sub create_read_file {
	# 引数を取得
	my ($file, $seed_id) = @_;
	
	# ファイルハンドルが定義されている場合は0を返す
	return(0) if defined($file);
	
	# シードIDを更新
	$seed_id++;
	$_[1] = $seed_id;
	
	# リード配列ファイルを作成
	open($_[0]->[0], ">", "$opt{o}/seed${seed_id}_upstream_r1.fa") or &exception::error("failed to make file: $opt{o}/seed${seed_id}_upstream_r1.fa");
	open($_[0]->[1], ">", "$opt{o}/seed${seed_id}_upstream_r2.fa") or &exception::error("failed to make file: $opt{o}/seed${seed_id}_upstream_r2.fa");
	open($_[0]->[2], ">", "$opt{o}/seed${seed_id}_downstream_r1.fa") or &exception::error("failed to make file: $opt{o}/seed${seed_id}_downstream_r1.fa");
	open($_[0]->[3], ">", "$opt{o}/seed${seed_id}_downstream_r2.fa") or &exception::error("failed to make file: $opt{o}/seed${seed_id}_downstream_r2.fa");
	return(1);
}

# サブルーチンを追加
### 編集範囲 終了 ###

## ここから例外処理のパッケージ ##
package exception;

# ヘルプ表示
sub help {
	print STDERR "$software ";
	print STDERR $specified_command ? $specified_command : $version;
	print STDERR "\n\nFunctions:\n  $note\n\nUsage:\n  $software ";
	if (!$specified_command and @command_list) {
		print STDERR "<command>\n";
		print STDERR "\nCommand:\n";
		foreach (@command_list) {print STDERR "  $_\t$command{$_}\n";}
	}
	else {
		print STDERR "$specified_command " if $specified_command;
		print STDERR "[options] " if @option_list;
		print STDERR "$usage\n";
		print STDERR "\nOptions:\n" if @option_list;
		foreach (@option_list) {print STDERR "  -$_\t$option{$_}\n";}
	}
	exit(0);
}

# エラー表示
sub error {
	print STDERR $software;
	print STDERR " $specified_command" if $specified_command;
	print STDERR ": Error: $_[0]";
	print STDERR ": $_[1] line $." if $_[1];
	print STDERR "\n";
	#threads->tid or map {$_->detach} threads->list;	# threadsモジュールを使用する場合はアンコメント
	exit(1);
}

# 注意表示
sub caution {
	print STDERR $software;
	print STDERR " $specified_command" if $specified_command;
	print STDERR ": Caution: $_[0]";
	print STDERR ": $_[1] line $." if $_[1];
	print STDERR "\n";
	return(1);
}

### 編集範囲 開始 ###
# サブルーチンを追加

# パッケージを追加
### 編集範囲 終了 ###
