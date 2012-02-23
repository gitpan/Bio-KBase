use 5.006;
use strict;
use warnings;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME                => 'Bio::KBase',
    AUTHOR              => q{Robert Olson <olson@mcs.anl.gov>},
    VERSION_FROM        => 'lib/Bio/KBase.pm',
    ABSTRACT_FROM       => 'lib/Bio/KBase.pm',
    ($ExtUtils::MakeMaker::VERSION >= 6.3002
      ? ('LICENSE'=> 'perl')
      : ()),
    PL_FILES            => {},
    PREREQ_PM => {
        'Test::More' => 0,
	'JSON::RPC::Client' => 0,
    },
    BUILD_REQUIRES => {
	'String::Random' => 0,
    },
    dist                => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', },
    clean               => { FILES => 'Bio-KBase-*' },
);
