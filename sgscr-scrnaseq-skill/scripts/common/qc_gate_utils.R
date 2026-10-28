# qc_gate_utils.R — Quality_Gate 三态评估通用函数
# Validates: Requirements 8.1-8.3 (P7)
#
# Feature: sgscr-scrnaseq-skill, Property 7: 对任意评估输入，给出且仅给出三态之一：
#   PASS / FAIL(附原因) / NEEDS_HUMAN_REVIEW(附待复查产物路径)，三态互斥且完整。
#
# 判定优先级（保证互斥完整）：
#   1. auto_pass == FALSE            → FAIL(detail)
#   2. auto_pass == TRUE 且 human_review_file 非空 → NEEDS_HUMAN_REVIEW(review_file)
#   3. 其余                          → PASS(detail)

sgs_qc_gate <- function(auto_pass, detail = "", human_review_file = NULL) {
  if (!isTRUE(auto_pass)) {
    list(status = "FAIL",
         detail = if (nzchar(detail)) detail else "auto gate failed (no detail given)")
  } else if (!is.null(human_review_file) && nzchar(human_review_file)) {
    list(status = "NEEDS_HUMAN_REVIEW",
         detail = sprintf("auto checks passed; human review required: %s; %s",
                          human_review_file, detail))
  } else {
    list(status = "PASS", detail = detail)
  }
}
