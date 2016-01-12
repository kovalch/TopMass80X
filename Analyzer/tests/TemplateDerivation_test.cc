/*
 * TemplateDerivation_test.cc
 *
 *  Created on: 10.12.2015
 *      Author: stadie
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE TemplateDerivationTest
#include <boost/test/unit_test.hpp>
#include "TemplateDerivation.h"

BOOST_AUTO_TEST_CASE(template_derivation_test_TemplateName) {
  BOOST_CHECK_EQUAL(TemplateDerivation::templateName(172.5, 1.02),
                    "jes102mass1725");
  BOOST_CHECK_EQUAL(TemplateDerivation::templateName(168.0, 0.98),
                    "jes098mass1680");
  BOOST_CHECK_EQUAL(
      TemplateDerivation::templateName(166.5, 0.95999999),
      "jes096mass1665");
}
