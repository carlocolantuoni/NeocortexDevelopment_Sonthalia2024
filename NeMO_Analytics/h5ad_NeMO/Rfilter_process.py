import ast
import re
import sys

class RFilterTranslator(ast.NodeVisitor):
    def __init__(self):
        self.result = []
        self.python_version = sys.version_info
    
    def visit(self, node):
        if node is None:
            return "None"
        result = super().visit(node)
        return result
    
    def visit_Constant(self, node):
        return repr(node.value)
    
    def visit_Str(self, node):
        return repr(node.s)
    
    def visit_Num(self, node):
        return str(node.n)
    
    def visit_NameConstant(self, node):
        return str(node.value)
        
    def visit_Name(self, node):
        column_name = node.id
        if self._needs_backticks(column_name):
            return f"`{column_name}`"
        return column_name
    
    def _needs_backticks(self, name):
        special_chars = ['.', ' ', '-', '+', '*', '/', '(', ')', '[', ']']
        return any(char in name for char in special_chars)
    
    def visit_Attribute(self, node):
        obj = self.visit(node.value)
        attr = node.attr
        return f"`{obj}.{attr}`"

    def visit_Call(self, node):
        if isinstance(node.func, ast.Name) and node.func.id == 'list':
            elements = []
            for arg in node.args:
                element_value = self.visit(arg)
                elements.append(element_value)
            return f"[{', '.join(elements)}]"
        elif isinstance(node.func, ast.Attribute):
            obj = self.visit(node.func.value)
            method = node.func.attr
            args = []
            for arg in node.args:
                arg_value = self.visit(arg)
                args.append(str(arg_value))
            if args:
                return f"{obj}.{method}({', '.join(args)})"
            else:
                return f"{obj}.{method}()"
        return self.generic_visit(node)

    def visit_Compare(self, node):
        left = self.visit(node.left)
        
        if len(node.ops) == 1 and len(node.comparators) == 1:
            op = node.ops[0]
            right = self.visit(node.comparators[0])
            if isinstance(op, ast.Eq):
                return f"({left} == {right})"
            elif isinstance(op, ast.NotEq):
                return f"({left} != {right})"
            elif isinstance(op, ast.Lt):
                return f"({left} < {right})"
            elif isinstance(op, ast.LtE):
                return f"({left} <= {right})"
            elif isinstance(op, ast.Gt):
                return f"({left} > {right})"
            elif isinstance(op, ast.GtE):
                return f"({left} >= {right})"
        return "UNKNOWN_COMPARE"
    
    def visit_BoolOp(self, node):
        if isinstance(node.op, ast.And):
            op_str = " and "
        elif isinstance(node.op, ast.Or):
            op_str = " or "
        else:
            op_str = " unknown_op "
        values = [self.visit(value) for value in node.values]
        return f"({op_str.join(values)})"
    
    def visit_UnaryOp(self, node):
        if isinstance(node.op, ast.Not):
            operand = self.visit(node.operand)
            return f"(not {operand})"
        return self.generic_visit(node)
    
    def visit_List(self, node):
        elements = []
        for elt in node.elts:
            element_value = self.visit(elt)
            elements.append(str(element_value))
        return f"[{', '.join(elements)}]"
    
    def generic_visit(self, node):
        return f"UNHANDLED_{type(node).__name__}"

def find_c_function_end(text, start_pos):
    i = start_pos
    paren_count = 0
    in_quote = False
    quote_char = None
    while i < len(text):
        char = text[i]
        if not in_quote:
            if char in ["'", '"']:
                in_quote = True
                quote_char = char
            elif char == '(':
                paren_count += 1
            elif char == ')':
                paren_count -= 1
                if paren_count == 0:
                    return i  # 找到了真正的结束位置
        else:
            if char == quote_char:
                # 检查是否是转义的引号
                if i > 0 and text[i-1] != '\\':
                    in_quote = False
                    quote_char = None
        i += 1
    return -1

def replace_c_functions(text):
    new_text = ""
    i = 0
    while i < len(text):
        match = re.search(r'(\w+)\s*%in%\s*c\s*\(', text[i:])
        if not match:
            new_text += text[i:]
            break
        match_start = i + match.start()
        new_text += text[i:match_start]
        var_name = match.group(1)
        c_start = i + match.end() - 1  
        c_end = find_c_function_end(text, c_start)
        
        if c_end != -1:
            c_content = text[c_start + 1:c_end]
            new_text += f"{var_name}.isin([{c_content}])"
            i = c_end + 1
        else:
            new_text += text[match_start:match_start + len(match.group(0))]
            i = match_start + len(match.group(0))
    
    return new_text

def preprocess_r_to_python(expression):
    expression = replace_c_functions(expression)
    expression = re.sub(
        r"grepl\s*\(\s*['\"]([^'\"]*)['\"],\s*(\w+),\s*ignore\.case\s*=\s*TRUE\s*\)",
        r"\2.str.contains('\1', case=False)",
        expression
    )
    expression = re.sub(
        r"grepl\s*\(\s*['\"]([^'\"]*)['\"],\s*(\w+)\s*\)",
        r"\2.str.contains('\1')",
        expression
    )
    expression = re.sub(r'mean\s*\(\s*(\w+)\s*\)', r'\1.mean()', expression)
    expression = re.sub(r'is\.na\s*\(\s*(\w+)\s*\)', r'\1.isna()', expression)
    expression = re.sub(r'(\w+)\s*%in%\s*list\((.*?)\)', r'\1.isin([\2])', expression)
    expression = expression.replace('&', ' and ')
    expression = expression.replace('|', ' or ')
    expression = expression.replace('!', ' not ')
    return expression.strip()


def translate_r_filter(r_expression):
    r_expression = preprocess_r_to_python(r_expression)
    try:
        tree = ast.parse(r_expression, mode='eval')
        translator = RFilterTranslator()
        result = translator.visit(tree.body)
        return result
    except SyntaxError as e:
        raise ValueError(f"Unable to parse R expression:: {r_expression}\nError: {e}")
        return None

if __name__ == "__main__":
    complex_tests = [
        "cell_type %in% c('Neuron', 'Astrocyte') & nCount_RNA > 1000 & !is.na(sample_id)",
        "nFeature_RNA >= 200 & nFeature_RNA <= 2500 & percent.mt < 5",
        "treatment %in% c('drugA', 'drugB') & response == 'responder'",
        "!cell_type %in% c('Doublet(maybe)', 'Unknown') & is.na(batch)",
        "grepl('Tumor', tissue) & !grepl('Metastasis', tissue)",
        "!(is.na(sample_id) | sample_id == '') & percent.ribo < 10",
        "geneX_expr > mean(geneX_expr) & cluster %in% c('1', '2', '3')"
    ]

    simple_tests = [
        "cell_type %in% c('B cell', 'T cell', 'NK cell')",
        "!cell_type %in% c('Doublet', 'Unknown')", 
        "nCount_RNA > 1000 & nFeature_RNA < 6000 & percent.mt < 10",
        "disease_subtype %in% c('HCC', 'ICC') | treatment_status == 'pre-treatment'"
    ]

    for test in complex_tests:
        try:
            result = translate_r_filter(test)
            print(f"Success: {test} → {result}")
        except Exception as e:
            print(f"Failed: {test} → {e}")