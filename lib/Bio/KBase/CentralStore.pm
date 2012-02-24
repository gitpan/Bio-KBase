package Bio::KBase::CentralStore;

use JSON::RPC::Client;
use strict;
use Data::Dumper;
use URI;

=head1 NAME

Bio::KBase::CentralStore

=head1 DESCRIPTION



=cut

=head2 $store = Bio::KBase::CentralStore->new($url)

Create a central store client.

=cut

sub new
{
    my($class, $url) = @_;

    $url ||= "http://bio-data-1.mcs.anl.gov/services/cdmi_api";
    
    my $self = {
	client => JSON::RPC::Client->new,
	url => $url,
    };
    return bless $self, $class;
}




=head2 $result = fids_to_annotations(fids)

This routine takes as input a list of fids.  It retrieves the existing
annotations for each fid, including the text of the annotation, who
made the annotation and when (as seconds from the epoch).

=cut

sub fids_to_annotations
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.fids_to_annotations",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking fids_to_annotations: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking fids_to_annotations: " . $self->{client}->status_line;
    }
}



=head2 $result = fids_to_functions(fids)

This routine takes as input a list of fids and returns a mapping
from the fids to their assigned functions.

=cut

sub fids_to_functions
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.fids_to_functions",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking fids_to_functions: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking fids_to_functions: " . $self->{client}->status_line;
    }
}



=head2 $result = fids_to_literature(fids)

We try to associate features and publications, when the publications constitute
supporting evidence of the function.  We connect a paper to a feature when
we believe that an "expert" has asserted that the function of the feature
is basically what we have associated with the feature.  Thus, we might
attach a paper reporting the crystal structure of a protein, even though
the paper is clearly not the paper responsible for the original characterization.
Our position in this matter is somewhat controversial, but we are seeking to
characterize some assertions as relatively solid, and this strategy seems to
support that goal.  Please note that we certainly wish we could also
capture original publications, and when experts can provide those
connections, we hope that they will help record the associations.

=cut

sub fids_to_literature
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.fids_to_literature",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking fids_to_literature: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking fids_to_literature: " . $self->{client}->status_line;
    }
}



=head2 $result = fids_to_protein_families(fids)

Kbase supports the creation and maintence of protein families.  Each family is intended to contain a set
of isofunctional homologs.  Currently, the families are collections of translations
of features, rather than of just protein sequences (represented by md5s, for example).
fids_to_protein_families supports access to the features that have been grouped into a family.
Ideally, each feature in a family would have the same assigned function.  This is not
always true, but probably should be.

=cut

sub fids_to_protein_families
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.fids_to_protein_families",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking fids_to_protein_families: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking fids_to_protein_families: " . $self->{client}->status_line;
    }
}



=head2 $result = fids_to_roles(fids)

Given a feature, one can get the set of roles it implements using fid_to_roles.
Remember, a protein can be multifunctional -- implementing several roles.
This can occur due to fusions or to broad specificity of substrate.

=cut

sub fids_to_roles
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.fids_to_roles",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking fids_to_roles: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking fids_to_roles: " . $self->{client}->status_line;
    }
}



=head2 $result = fids_to_subsystems(fids)

fids in subsystems normally have somewhat more reliable assigned functions than
those not in subsystems.  Hence, it is common to ask "Is this protein-encoding gene
included in any subsystems?"   fids_to_subsystems can be used to see which subsystems
contain a fid (or, you can submit as input a set of fids and get the subsystems for each).

=cut

sub fids_to_subsystems
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.fids_to_subsystems",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking fids_to_subsystems: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking fids_to_subsystems: " . $self->{client}->status_line;
    }
}



=head2 $result = fids_to_co_occurring_fids(fids)

One of the most powerful clues to function relates to conserved clusters of genes on
the chromosome (in prokaryotic genomes).  We have attempted to record pairs of genes
that tend to occur close to one another on the chromosome.  To meaningfully do this,
we need to construct similarity-based mappings between genes in distinct genomes.
We have constructed such mappings for many (but not all) genomes maintained in the
Kbase CS.  The prokaryotic geneomes in the CS are grouped into OTUs by ribosomal
RNA (genomes within a single OTU have SSU rRNA that is greater than 97% identical).
If two genes occur close to one another (i.e., corresponding genes occur close
to one another), then we assign a score, which is the number of distinct OTUs
in which such clustering is detected.  This allows one to normalize for situations
in which hundreds of corresponding genes are detected, but they all come from
very closely related genomes.

The significance of the score relates to the number of genomes in the database.
We recommend that you take the time to look at a set of scored pairs and determine
approximately what percentage appear to be actually related for a few cutoff values.

=cut

sub fids_to_co_occurring_fids
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.fids_to_co_occurring_fids",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking fids_to_co_occurring_fids: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking fids_to_co_occurring_fids: " . $self->{client}->status_line;
    }
}



=head2 $result = fids_to_locations(fids)

A "location" is a sequence of "regions".  A region is a contiguous set of bases
in a contig.  We work with locations in both the string form and as structures.
fids_to_locations takes as input a list of fids.  For each fid, a structured location
is returned.  The location is a list of regions; a region is given as a pointer to
a list containing

             the contig,
             the beginning base in the contig (from 1).
             the strand (+ or -), and
             the length

Note that specifying a region using these 4 values allows you to represent a single
base-pair region on either strand unambiguously (which giving begin/end pairs does
not achieve).

=cut

sub fids_to_locations
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.fids_to_locations",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking fids_to_locations: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking fids_to_locations: " . $self->{client}->status_line;
    }
}



=head2 $result = locations_to_fids(region_of_dna_strings)

It is frequently the case that one wishes to look up the genes that
occur in a given region of a contig.  Location_to_fids can be used to extract
such sets of genes for each region in the input set of regions.  We define a gene
as "occuring" in a region if the location of the gene overlaps the designated region.

=cut

sub locations_to_fids
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.locations_to_fids",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking locations_to_fids: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking locations_to_fids: " . $self->{client}->status_line;
    }
}



=head2 $result = locations_to_dna_sequences(locations)

locations_to_dna_sequences takes as input a list of locations (each in the form of
a list of regions).  The routine constructs 2-tuples composed of

     [the input location,the dna string]

The returned DNA string is formed by concatenating the DNA for each of the
regions that make up the location.

=cut

sub locations_to_dna_sequences
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.locations_to_dna_sequences",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking locations_to_dna_sequences: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking locations_to_dna_sequences: " . $self->{client}->status_line;
    }
}



=head2 $result = proteins_to_fids(proteins)

proteins_to_fids takes as input a list of proteins (i.e., a list of md5s) and
returns for each a set of protein-encoding fids that have the designated
sequence as their translation.  That is, for each sequence, the returned fids will
be the entire set (within Kbase) that have the sequence as a translation.

=cut

sub proteins_to_fids
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.proteins_to_fids",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking proteins_to_fids: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking proteins_to_fids: " . $self->{client}->status_line;
    }
}



=head2 $result = proteins_to_protein_families(proteins)

Protein families contain a set of isofunctional homologs.  proteins_to_protein_families
can be used to look up is used to get the set of protein_families containing a specified protein.
For performance reasons, you can submit a batch of proteins (i.e., a list of proteins),
and for each input protein, you get back a set (possibly empty) of protein_families.
Specific collections of families (e.g., FIGfams) usually require that a protein be in
at most one family.  However, we will be integrating protein families from a number of
sources, and so a protein can be in multiple families.

=cut

sub proteins_to_protein_families
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.proteins_to_protein_families",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking proteins_to_protein_families: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking proteins_to_protein_families: " . $self->{client}->status_line;
    }
}



=head2 $result = proteins_to_literature(proteins)

The routine proteins_to_literature can be used to extract the list of papers
we have associated with specific protein sequences.  The user should note that
in many cases the association of a paper with a protein sequence is not precise.
That is, the paper may actually describe a closely-related protein (that may
not yet even be in a sequenced genome).  Annotators attempt to use best
judgement when associating literature and proteins.  Publication references
include [pubmed ID,URL for the paper, title of the paper].  In some cases,
the URL and title are omitted.  In theory, we can extract them from PubMed
and we will attempt to do so.

=cut

sub proteins_to_literature
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.proteins_to_literature",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking proteins_to_literature: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking proteins_to_literature: " . $self->{client}->status_line;
    }
}



=head2 $result = proteins_to_functions(proteins)

The routine proteins_to_functions allows users to access functions associated with
specific protein sequences.  The input proteins are given as a list of MD5 values
(these MD5 values each correspond to a specific protein sequence).  For each input
MD5 value, a list of [feature-id,function] pairs is constructed and returned.
Note that there are many cases in which a single protein sequence corresponds
to the translation associated with multiple protein-encoding genes, and each may
have distinct functions (an undesirable situation, we grant).

This function allows you to access all of the functions assigned (by all annotation
groups represented in Kbase) to each of a set of sequences.

=cut

sub proteins_to_functions
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.proteins_to_functions",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking proteins_to_functions: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking proteins_to_functions: " . $self->{client}->status_line;
    }
}



=head2 $result = proteins_to_roles(proteins)

The routine proteins_to_roles allows a user to gather the set of functional
roles that are associated with specifc protein sequences.  A single protein
sequence (designated by an MD5 value) may have numerous associated functions,
since functions are treated as an attribute of the feature, and multiple
features may have precisely the same translation.  In our experience,
it is not uncommon, even for the best annotation teams, to assign
distinct functions (and, hence, functional roles) to identical
protein sequences.

For each input MD5 value, this routine gathers the set of features (fids)
that share the same sequence, collects the associated functions, expands
these into functional roles (for multi-functional proteins), and returns
the set of roles that results.

Note that, if the user wishes to see the specific features that have the
assigned fiunctional roles, they should use proteins_to_functions instead (it
returns the fids associated with each assigned function).

=cut

sub proteins_to_roles
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.proteins_to_roles",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking proteins_to_roles: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking proteins_to_roles: " . $self->{client}->status_line;
    }
}



=head2 $result = roles_to_proteins(roles)

roles_to_proteins can be used to extract the set of proteins (designated by MD5 values)
that currently are believed to implement a given role.  Note that the proteins
may be multifunctional, meaning that they may be implementing other roles, as well.

=cut

sub roles_to_proteins
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.roles_to_proteins",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking roles_to_proteins: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking roles_to_proteins: " . $self->{client}->status_line;
    }
}



=head2 $result = roles_to_subsystems(roles)

roles_to_subsystems can be used to access the set of subsystems that include
specific roles. The input is a list of roles (i.e., role descriptions), and a mapping
is returned as a hash with key role description and values composed of sets of susbsystem names.

=cut

sub roles_to_subsystems
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.roles_to_subsystems",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking roles_to_subsystems: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking roles_to_subsystems: " . $self->{client}->status_line;
    }
}



=head2 $result = roles_to_protein_families(roles)

roles_to_protein_families can be used to locate the protein families containing
features that have assigned functions implying that they implement designated roles.
Note that for any input role (given as a role description), you may have a set
of distinct protein_families returned.

=cut

sub roles_to_protein_families
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.roles_to_protein_families",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking roles_to_protein_families: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking roles_to_protein_families: " . $self->{client}->status_line;
    }
}



=head2 $result = fids_to_coexpressed_fids(fids)

The routine fids_to_coexpressed_fids returns (for each input fid) a
list of features that appear to be coexpressed.  That is,
for an input fid, we determine the set of fids from the same genome that
have Pearson Correlation Coefficients (based on normalized expression data)
greater than 0.5 or less than -0.5.

=cut

sub fids_to_coexpressed_fids
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.fids_to_coexpressed_fids",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking fids_to_coexpressed_fids: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking fids_to_coexpressed_fids: " . $self->{client}->status_line;
    }
}



=head2 $result = protein_families_to_fids(protein_families)

protein_families_to_fids can be used to access the set of fids represented by each of
a set of protein_families.  We define protein_families as sets of fids (rather than sets
of MD5s.  This may, or may not, be a mistake.

=cut

sub protein_families_to_fids
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.protein_families_to_fids",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking protein_families_to_fids: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking protein_families_to_fids: " . $self->{client}->status_line;
    }
}



=head2 $result = protein_families_to_proteins(protein_families)

protein_families_to_proteins can be used to access the set of proteins (i.e., the set of MD5 values)
represented by each of a set of protein_families.  We define protein_families as sets of fids (rather than sets
           of MD5s.  This may, or may not, be a mistake.

=cut

sub protein_families_to_proteins
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.protein_families_to_proteins",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking protein_families_to_proteins: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking protein_families_to_proteins: " . $self->{client}->status_line;
    }
}



=head2 $result = protein_families_to_functions(protein_families)

protein_families_to_functions can be used to extract the set of functions assigned to the fids
that make up the family.  Each input protein_family is mapped to a family function.

=cut

sub protein_families_to_functions
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.protein_families_to_functions",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking protein_families_to_functions: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking protein_families_to_functions: " . $self->{client}->status_line;
    }
}



=head2 $result = protein_families_to_co_occurring_families(protein_families)

Since we accumulate data relating to the co-occurrence (i.e., chromosomal
clustering) of genes in prokaryotic genomes,  we can note which pairs of genes tend to co-occur.
From this data, one can compute the protein families that tend to co-occur (i.e., tend to
cluster on the chromosome).  This allows one to formulate conjectures for unclustered pairs, based
on clustered pairs from the same protein_families.

=cut

sub protein_families_to_co_occurring_families
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.protein_families_to_co_occurring_families",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking protein_families_to_co_occurring_families: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking protein_families_to_co_occurring_families: " . $self->{client}->status_line;
    }
}



=head2 $result = co_occurrence_evidence(pairs_of_fids)

co-occurence_evidence is used to retrieve the detailed pairs of genes that go into the
computation of co-occurence scores.  The scores reflect an estimate of the number of distinct OTUs that
contain an instance of a co-occuring pair.  This routine returns as evidence a list of all the pairs that
went into the computation.

The input to the computation is a list of pairs for which evidence is desired.

The returned output is a list of elements. one for each input pair.  Each output element
is a 2-tuple: the input pair and the evidence for the pair.  The evidence is a list of pairs of
fids that are believed to correspond to the input pair.

=cut

sub co_occurrence_evidence
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.co_occurrence_evidence",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking co_occurrence_evidence: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking co_occurrence_evidence: " . $self->{client}->status_line;
    }
}



=head2 $result = contigs_to_sequences(contigs)

contigs_to_sequences is used to access the DNA sequence associated with each of a set
of input contigs.  It takes as input a set of contig IDs (from which the genome can be determined) and
produces a mapping from the input IDs to the returned DNA sequence in each case.

=cut

sub contigs_to_sequences
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.contigs_to_sequences",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking contigs_to_sequences: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking contigs_to_sequences: " . $self->{client}->status_line;
    }
}



=head2 $result = contigs_to_lengths(contigs)

In some cases, one wishes to know just the lengths of the contigs, rather than their
actual DNA sequence (e.g., suppose that you wished to know if a gene boundary occured within
100 bp of the end of the contig).  To avoid requiring a user to access the entire DNA sequence,
we offer the ability to retrieve just the contig lengths.  Input to the routine is a list of contig IDs.
The routine returns a mapping from contig IDs to lengths

=cut

sub contigs_to_lengths
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.contigs_to_lengths",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking contigs_to_lengths: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking contigs_to_lengths: " . $self->{client}->status_line;
    }
}



=head2 $result = contigs_to_md5s(contigs)

contigs_to_md5s can be used to acquire MD5 values for each of a list of contigs.
The quickest way to determine whether two contigs are identical is to compare their
associated MD5 values, eliminating the need to retrieve the sequence of each and compare them.

The routine takes as input a list of contig IDs.  The output is a mapping
from contig ID to MD5 value.

=cut

sub contigs_to_md5s
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.contigs_to_md5s",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking contigs_to_md5s: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking contigs_to_md5s: " . $self->{client}->status_line;
    }
}



=head2 $result = md5s_to_genomes(md5s)

md5s to genomes is used to get the genomes associated with each of a list of input md5 values.

           The routine takes as input a list of MD5 values.  It constructs a mapping from each input
           MD5 value to a list of genomes that share the same MD5 value.

           The MD5 value for a genome is independent of the names of contigs and the case of the DNA sequence
           data.

=cut

sub md5s_to_genomes
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.md5s_to_genomes",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking md5s_to_genomes: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking md5s_to_genomes: " . $self->{client}->status_line;
    }
}



=head2 $result = genomes_to_md5s(genomes)

The routine genomes_to_md5s can be used to look up the MD5 value associated with each of
a set of genomes.  The MD5 values are computed when the genome is loaded, so this routine
just retrieves the precomputed values.

Note that the MD5 value of a genome is independent of the contig names and case of the
DNA sequences that make up the genome.

=cut

sub genomes_to_md5s
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.genomes_to_md5s",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking genomes_to_md5s: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking genomes_to_md5s: " . $self->{client}->status_line;
    }
}



=head2 $result = genomes_to_contigs(genomes)

The routine genomes_to_contigs can be used to retrieve the IDs of the contigs
associated with each of a list of input genomes.  The routine constructs a mapping
from genome ID to the list of contigs included in the genome.

=cut

sub genomes_to_contigs
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.genomes_to_contigs",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking genomes_to_contigs: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking genomes_to_contigs: " . $self->{client}->status_line;
    }
}



=head2 $result = genomes_to_fids(genomes, types_of_fids)

genomes_to_fids is used to get the fids included in specific genomes.  It
is often the case that you want just one or two types of fids -- hence, the
types_of_fids argument.

=cut

sub genomes_to_fids
{
    my($self, @args) = @_;

    @args == 2 or die "Invalid argument count (expecting 2)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.genomes_to_fids",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking genomes_to_fids: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking genomes_to_fids: " . $self->{client}->status_line;
    }
}



=head2 $result = genomes_to_taxonomies(genomes)

The routine genomes_to_taxonomies can be used to retrieve taxonomic information for
each of a list of input genomes.  For each genome in the input list of genomes, a list of
taxonomic groups is returned.  Kbase will use the groups maintained by NCBI.  For an NCBI
taxonomic string like

     cellular organisms;
     Bacteria;
     Proteobacteria;
     Gammaproteobacteria;
     Enterobacteriales;
     Enterobacteriaceae;
     Escherichia;
     Escherichia coli

associated with the strain 'Escherichia coli 1412', this routine would return a list of these
taxonomic groups:


     ['Bacteria',
      'Proteobacteria',
      'Gammaproteobacteria',
      'Enterobacteriales',
      'Enterobacteriaceae',
      'Escherichia',
      'Escherichia coli',
      'Escherichia coli 1412'
     ]

That is, the initial "cellular organisms" has been deleted, and the strain ID has
been added as the last "grouping".

The output is a mapping from genome IDs to lists of the form shown above.

=cut

sub genomes_to_taxonomies
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.genomes_to_taxonomies",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking genomes_to_taxonomies: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking genomes_to_taxonomies: " . $self->{client}->status_line;
    }
}



=head2 $result = genomes_to_subsystems(genomes)

A user can invoke genomes_to_subsystems to rerieve the names of the subsystems
relevant to each genome.  The input is a list of genomes.  The output is a mapping
from genome to a list of 2-tuples, where each 2-tuple give a variant code and a
subsystem name.  Variant codes of -1 (or *-1) amount to assertions that the
genome contains no active variant.  A variant code of 0 means "work in progress",
and presence or absence of the subsystem in the genome should be undetermined.

=cut

sub genomes_to_subsystems
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.genomes_to_subsystems",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking genomes_to_subsystems: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking genomes_to_subsystems: " . $self->{client}->status_line;
    }
}



=head2 $result = subsystems_to_genomes(subsystems)

The routine subsystems_to_genomes is used to determine which genomes are in
specified subsystems.  The input is the list of subsystem names of interest.
The output is a map from the subsystem names to lists of 2-tuples, where each 2-tuple is
a [variant-code,genome ID] pair.

=cut

sub subsystems_to_genomes
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.subsystems_to_genomes",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking subsystems_to_genomes: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking subsystems_to_genomes: " . $self->{client}->status_line;
    }
}



=head2 $result = subsystems_to_fids(subsystems, genomes)

The routine subsystems_to_fids allows the user to map subsystem names into the fids that
occur in genomes in the subsystems.  Specifically, the input is a list of subsystem names.
What is returned is a mapping from subsystem names to a "genome-mapping".  The genome-mapping
takes genome IDs to 2-tuples that capture the variant code of the genome and the fids from
the genome that are included in the subsystem.

=cut

sub subsystems_to_fids
{
    my($self, @args) = @_;

    @args == 2 or die "Invalid argument count (expecting 2)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.subsystems_to_fids",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking subsystems_to_fids: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking subsystems_to_fids: " . $self->{client}->status_line;
    }
}



=head2 $result = subsystems_to_roles(subsystems, aux)

The routine subsystem_to_roles is used to determine the role descriptions that
occur in a subsystem.  The input is a list of subsystem names.  A map is returned connecting
subsystem names to lists of roles.  'aux' is a boolean variable.  If it is 0, auxiliary roles
are not returned.  If it is 1, they are returned.

=cut

sub subsystems_to_roles
{
    my($self, @args) = @_;

    @args == 2 or die "Invalid argument count (expecting 2)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.subsystems_to_roles",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking subsystems_to_roles: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking subsystems_to_roles: " . $self->{client}->status_line;
    }
}



=head2 $result = subsystems_to_spreadsheets(subsystems, genomes)

The subsystem_to_spreadsheet routine allows a user to extract the subsystem spreadsheets for
a specified set of subsystem names.  In the returned output, each subsystem is mapped
to a hash that takes as input a genome ID and maps it to the "row" for the genome in the subsystem.
The "row" is itself a 2-tuple composed of the variant code, and a mapping from role descriptions to
lists of fids.  We suggest writing a simple test script to get, say, the subsystem named
'Histidine Degradation', extracting the spreadsheet, and then using something like Dumper to make
sure that it all makes sense.

=cut

sub subsystems_to_spreadsheets
{
    my($self, @args) = @_;

    @args == 2 or die "Invalid argument count (expecting 2)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.subsystems_to_spreadsheets",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking subsystems_to_spreadsheets: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking subsystems_to_spreadsheets: " . $self->{client}->status_line;
    }
}



=head2 $result = all_roles_used_in_models()

The all_roles_used_in_models allows a user to access the set of roles that are included in current models.  This is
important.  There are far fewer roles used in models than overall.  Hence, the returned set represents
the minimal set we need to clean up in order to properly support modeling.

=cut

sub all_roles_used_in_models
{
    my($self, @args) = @_;

    @args == 0 or die "Invalid argument count (expecting 0)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.all_roles_used_in_models",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking all_roles_used_in_models: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking all_roles_used_in_models: " . $self->{client}->status_line;
    }
}



=head2 $result = complexes_to_complex_data(complexes)



=cut

sub complexes_to_complex_data
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.complexes_to_complex_data",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking complexes_to_complex_data: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking complexes_to_complex_data: " . $self->{client}->status_line;
    }
}



=head2 $result = genomes_to_genome_data(genomes)



=cut

sub genomes_to_genome_data
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.genomes_to_genome_data",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking genomes_to_genome_data: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking genomes_to_genome_data: " . $self->{client}->status_line;
    }
}



=head2 $result = fids_to_regulon_data(fids)



=cut

sub fids_to_regulon_data
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.fids_to_regulon_data",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking fids_to_regulon_data: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking fids_to_regulon_data: " . $self->{client}->status_line;
    }
}



=head2 $result = fids_to_feature_data(fids)



=cut

sub fids_to_feature_data
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.fids_to_feature_data",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking fids_to_feature_data: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking fids_to_feature_data: " . $self->{client}->status_line;
    }
}



=head2 $result = equiv_sequence_assertions(proteins)

Different groups have made assertions of function for numerous protein sequences.
The equiv_sequence_assertions allows the user to gather function assertions from
all of the sources.  Each assertion includes a field indicating whether the person making
the assertion viewed themself as an "expert".  The routine gathers assertions for all
proteins having identical protein sequence.

=cut

sub equiv_sequence_assertions
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.equiv_sequence_assertions",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking equiv_sequence_assertions: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking equiv_sequence_assertions: " . $self->{client}->status_line;
    }
}



=head2 $result = fids_to_atomic_regulons(fids)

The fids_to_atomic_regulons allows one to map fids into regulons that contain the fids.
Normally a fid will be in at most one regulon, but we support multiple regulons.

=cut

sub fids_to_atomic_regulons
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.fids_to_atomic_regulons",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking fids_to_atomic_regulons: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking fids_to_atomic_regulons: " . $self->{client}->status_line;
    }
}



=head2 $result = atomic_regulons_to_fids(regulons)

The atomic_regulons_to_fids routine allows the user to access the set of fids that make up a regulon.
Regulons may arise from several sources; hence, fids can be in multiple regulons.

=cut

sub atomic_regulons_to_fids
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.atomic_regulons_to_fids",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking atomic_regulons_to_fids: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking atomic_regulons_to_fids: " . $self->{client}->status_line;
    }
}



=head2 $result = fids_to_protein_sequences(fids)

fids_to_protein_sequences allows the user to look up the amino acid sequences
corresponding to each of a set of fids.  You can also get the sequence from proteins (i.e., md5 values).
This routine saves you having to look up the md5 sequence and then accessing
the protein string in a separate call.

=cut

sub fids_to_protein_sequences
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.fids_to_protein_sequences",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking fids_to_protein_sequences: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking fids_to_protein_sequences: " . $self->{client}->status_line;
    }
}



=head2 $result = fids_to_proteins(fids)



=cut

sub fids_to_proteins
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.fids_to_proteins",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking fids_to_proteins: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking fids_to_proteins: " . $self->{client}->status_line;
    }
}



=head2 $result = fids_to_dna_sequences(fids)

fids_to_dna_sequences allows the user to look up the DNA sequences
corresponding to each of a set of fids.

=cut

sub fids_to_dna_sequences
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.fids_to_dna_sequences",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking fids_to_dna_sequences: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking fids_to_dna_sequences: " . $self->{client}->status_line;
    }
}



=head2 $result = roles_to_fids(roles, genomes)

A "function" is a set of "roles" (often called "functional roles");

                F1 / F2  (where F1 and F2 are roles)  is a function that implements
                          two functional roles in different domains of the protein.
                F1 @ F2 implements multiple roles through broad specificity
                F1; F2  is thought to implement F1 or f2 (uncertainty)

            You often wish to find the fids in one or more genomes that
            implement specific functional roles.  To do this, you can use
            roles_to_fids.

=cut

sub roles_to_fids
{
    my($self, @args) = @_;

    @args == 2 or die "Invalid argument count (expecting 2)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.roles_to_fids",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking roles_to_fids: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking roles_to_fids: " . $self->{client}->status_line;
    }
}



=head2 $result = reactions_to_complexes(reactions)

Reactions are thought of as being either spontaneous or implemented by
one or more Complexes.  Complexes connect to Roles.  Hence, the connection of fids
or roles to reactions goes through Complexes.

=cut

sub reactions_to_complexes
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.reactions_to_complexes",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking reactions_to_complexes: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking reactions_to_complexes: " . $self->{client}->status_line;
    }
}



=head2 $result = reaction_strings(reactions, name_parameter)

Reaction_strings are text strings that represent (albeit crudely)
the details of Reactions.

=cut

sub reaction_strings
{
    my($self, @args) = @_;

    @args == 2 or die "Invalid argument count (expecting 2)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.reaction_strings",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking reaction_strings: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking reaction_strings: " . $self->{client}->status_line;
    }
}



=head2 $result = roles_to_complexes(roles)

roles_to_complexes allows a user to connect Roles to Complexes,
from there, the connection exists to Reactions (although in the
actual ER-model model, the connection from Complex to Reaction goes through
ReactionComplex).  Since Roles also connect to fids, the connection between
fids and Reactions is induced.

The "name_parameter" can be 0, 1 or 'only'. If 1, then the compound name will 
be included with the ID in the output. If only, the compound name will be included 
instead of the ID. If 0, only the ID will be included. The default is 0.

=cut

sub roles_to_complexes
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.roles_to_complexes",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking roles_to_complexes: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking roles_to_complexes: " . $self->{client}->status_line;
    }
}



=head2 $result = fids_to_subsystem_data(fids)



=cut

sub fids_to_subsystem_data
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.fids_to_subsystem_data",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking fids_to_subsystem_data: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking fids_to_subsystem_data: " . $self->{client}->status_line;
    }
}



=head2 $result = representative(genomes)



=cut

sub representative
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.representative",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking representative: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking representative: " . $self->{client}->status_line;
    }
}



=head2 $result = otu_members(genomes)



=cut

sub otu_members
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.otu_members",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking otu_members: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking otu_members: " . $self->{client}->status_line;
    }
}



=head2 $result = fids_to_genomes(fids)



=cut

sub fids_to_genomes
{
    my($self, @args) = @_;

    @args == 1 or die "Invalid argument count (expecting 1)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.fids_to_genomes",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking fids_to_genomes: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking fids_to_genomes: " . $self->{client}->status_line;
    }
}



=head2 $result = text_search(input, start, count, entities)

text_search performs a search against a full-text index maintained 
for the CDMI. The parameter "input" is the text string to be searched for.
The parameter "entities" defines the entities to be searched. If the list
is empty, all indexed entities will be searched. The "start" and "count"
parameters limit the results to "count" hits starting at "start".

=cut

sub text_search
{
    my($self, @args) = @_;

    @args == 4 or die "Invalid argument count (expecting 4)";
    my $result = $self->{client}->call($self->{url}, {
	method => "CDMI_API.text_search",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    die "Error invoking text_search: " . $result->error_message;
	} else {
	    return @{$result->result};
	}
    } else {
	die "Error invoking text_search: " . $self->{client}->status_line;
    }
}




1;
