-- phpMyAdmin SQL Dump
-- version 4.0.4.1
-- http://www.phpmyadmin.net
--
-- Host: localhost
-- Generation Time: Mar 30, 2016 at 04:00 PM
-- Server version: 5.0.95
-- PHP Version: 5.3.3

SET SQL_MODE = "NO_AUTO_VALUE_ON_ZERO";
SET time_zone = "+00:00";


/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;

--
-- Database: `ENVDB`
--
CREATE DATABASE IF NOT EXISTS `ENVDB` DEFAULT CHARACTER SET latin1 COLLATE latin1_swedish_ci;
USE `ENVDB`;

-- --------------------------------------------------------

--
-- Table structure for table `feature_accessions`
--

CREATE TABLE IF NOT EXISTS `feature_accessions` (
  `feature_id` int(11) NOT NULL,
  `source` varchar(50) NOT NULL,
  `prefix` varchar(50) default NULL,
  `accession` varchar(100) NOT NULL,
  KEY `feature_id` (`feature_id`,`source`,`accession`),
  KEY `accession` (`accession`),
  KEY `source` (`source`),
  KEY `feature_id_2` (`feature_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `feature_annotations`
--

CREATE TABLE IF NOT EXISTS `feature_annotations` (
  `feature_id` int(11) unsigned NOT NULL default '0',
  `data_type_id` int(11) unsigned NOT NULL default '0',
  `value` varchar(256) NOT NULL default '',
  `edit_by` varchar(30) NOT NULL default 'CURRENT_USER',
  `date_edit` timestamp NOT NULL default CURRENT_TIMESTAMP,
  `rank` tinyint(4) default NULL,
  `source` varchar(100) default NULL,
  `is_current` tinyint(1) NOT NULL default '1',
  KEY `feature_id` (`feature_id`,`data_type_id`),
  KEY `value_data_type_id` (`value`,`data_type_id`),
  KEY `source` (`source`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `feature_evidence`
--

CREATE TABLE IF NOT EXISTS `feature_evidence` (
  `feature_id` int(9) NOT NULL,
  `feat_min` int(6) NOT NULL,
  `feat_max` int(6) NOT NULL,
  `program` varchar(100) NOT NULL,
  `ev_type` varchar(100) NOT NULL,
  `ev_accession` varchar(100) default NULL,
  `ev_min` int(6) default NULL,
  `ev_max` int(6) default NULL,
  `ev_length` int(6) default NULL,
  `score` varchar(255) NOT NULL,
  `ev_date` timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,
  KEY `feature_id` (`feature_id`,`ev_type`,`ev_accession`),
  KEY `program` (`program`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `import`
--

CREATE TABLE IF NOT EXISTS `import` (
  `feature_id` int(11) NOT NULL,
  `source` varchar(50) NOT NULL,
  `prefix` varchar(50) default NULL,
  `accession` varchar(100) NOT NULL,
  KEY `feature_id` (`feature_id`,`source`,`accession`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `proteomics_data`
--

CREATE TABLE IF NOT EXISTS `proteomics_data` (
  `feature_id` int(11) unsigned NOT NULL,
  `experiment` varchar(250) NOT NULL,
  `sample` varchar(250) NOT NULL,
  `filt_spectra` int(7) unsigned NOT NULL,
  `norm_value` float NOT NULL,
  KEY `feature_id` (`feature_id`,`experiment`),
  KEY `sample` (`sample`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `sequences`
--

CREATE TABLE IF NOT EXISTS `sequences` (
  `seq_id` int(11) unsigned NOT NULL auto_increment,
  `category` enum('read','contig','scaffold','closed') NOT NULL,
  `sequence` mediumtext NOT NULL,
  `date_inserted` date NOT NULL default '0000-00-00',
  `inserted_by` varchar(30) NOT NULL default '',
  `seq_length` mediumint(8) unsigned default NULL,
  `gc` float unsigned default NULL,
  `iscurrent` binary(1) NOT NULL default ' ',
  PRIMARY KEY  (`seq_id`),
  KEY `seq_length` (`seq_length`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Holds production sequence.' AUTO_INCREMENT=1 ;

-- --------------------------------------------------------

--
-- Table structure for table `sequence_accessions`
--

CREATE TABLE IF NOT EXISTS `sequence_accessions` (
  `seq_id` int(11) unsigned NOT NULL default '0',
  `seq_accession` varchar(50) NOT NULL default '',
  `seq_acc_source` varchar(50) NOT NULL default '',
  `seq_acc_description` tinytext,
  PRIMARY KEY  (`seq_id`,`seq_accession`),
  KEY `seq_acc_source` (`seq_acc_source`),
  KEY `seq_accession` (`seq_accession`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='Links public and private accessions and names to sequences';

-- --------------------------------------------------------

--
-- Table structure for table `sequence_annotations`
--

CREATE TABLE IF NOT EXISTS `sequence_annotations` (
  `seq_id` int(11) unsigned NOT NULL default '0',
  `data_type` varchar(50) NOT NULL,
  `value` varchar(100) NOT NULL default '',
  `edit_by` varchar(30) NOT NULL default '',
  `date_edit` timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,
  `rank` tinyint(4) unsigned NOT NULL default '0',
  KEY `sequence_id` (`seq_id`),
  KEY `data_type` (`data_type`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `sequence_features`
--

CREATE TABLE IF NOT EXISTS `sequence_features` (
  `feature_id` int(11) unsigned NOT NULL auto_increment,
  `feat_type_id` int(11) unsigned NOT NULL default '0',
  `is_current` tinyint(1) NOT NULL default '1',
  `feat_type` varchar(25) NOT NULL,
  `SO_term` mediumint(7) unsigned zerofill default NULL,
  `product` text,
  `inserted_by` varchar(30) NOT NULL default '',
  `date_inserted` date NOT NULL default '0000-00-00',
  PRIMARY KEY  (`feature_id`),
  KEY `feat_type_id` (`feat_type_id`),
  KEY `SO_term` (`SO_term`),
  KEY `feat_type` (`feat_type`),
  KEY `is_current` (`is_current`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 AUTO_INCREMENT=1 ;

-- --------------------------------------------------------

--
-- Table structure for table `sequence_mappings`
--

CREATE TABLE IF NOT EXISTS `sequence_mappings` (
  `seq_id` int(11) unsigned NOT NULL default '0',
  `seq_min` mediumint(8) unsigned NOT NULL default '0',
  `seq_max` mediumint(8) unsigned NOT NULL default '0',
  `subseq_id` int(11) unsigned NOT NULL default '0',
  `sub_min` mediumint(8) unsigned NOT NULL default '0',
  `sub_max` mediumint(8) unsigned NOT NULL default '0',
  `revcomp` binary(1) NOT NULL default '0',
  KEY `sequence_id` (`seq_id`,`subseq_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `sequence_sets`
--

CREATE TABLE IF NOT EXISTS `sequence_sets` (
  `set_id` mediumint(8) unsigned NOT NULL auto_increment,
  `name` varchar(50) NOT NULL default '',
  `description` varchar(255) NOT NULL default '',
  `ncbi_proj_id` varchar(15) default NULL,
  `ncbi_sample_id` varchar(15) default NULL,
  `WGS_accession` varchar(15) default NULL,
  `NCBI_TAXON_ID` int(11) default NULL,
  `is_current` bit(1) NOT NULL,
  PRIMARY KEY  (`set_id`),
  UNIQUE KEY `name` (`name`),
  KEY `ncbi_proj_id` (`ncbi_proj_id`),
  KEY `ncbi_sample_id` (`ncbi_sample_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 AUTO_INCREMENT=1 ;

-- --------------------------------------------------------

--
-- Table structure for table `seq_feat_mappings`
--

CREATE TABLE IF NOT EXISTS `seq_feat_mappings` (
  `seq_id` int(11) unsigned NOT NULL default '0',
  `feature_id` int(11) unsigned NOT NULL default '0',
  `feat_min` mediumint(8) NOT NULL default '0',
  `feat_max` mediumint(8) NOT NULL default '0',
  `strand` enum('1','0','-1') NOT NULL default '0',
  `phase` enum('0','1','2') default NULL,
  `min_partial` binary(1) NOT NULL default '0',
  `max_partial` binary(1) NOT NULL default '0',
  `pseudo` int(1) unsigned NOT NULL default '0' COMMENT '0 - not pseudo, 1 - single FS, 2 - single STOP, 3 - single indel, 4 - multiple, 5 - unspecified',
  KEY `sequence_id` (`seq_id`,`feature_id`),
  KEY `pseudo` (`pseudo`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `seq_set_link`
--

CREATE TABLE IF NOT EXISTS `seq_set_link` (
  `seq_id` mediumint(9) NOT NULL default '0',
  `set_id` mediumint(9) NOT NULL default '0',
  KEY `sequence_id` (`seq_id`,`set_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `transcriptomics_data`
--

CREATE TABLE IF NOT EXISTS `transcriptomics_data` (
  `feature_id` int(11) unsigned NOT NULL,
  `experiment` varchar(250) NOT NULL,
  `sample` varchar(250) NOT NULL,
  `abs_value` int(9) default NULL,
  `norm_value` float unsigned NOT NULL,
  KEY `feature_id` (`feature_id`,`experiment`),
  KEY `sample` (`sample`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
